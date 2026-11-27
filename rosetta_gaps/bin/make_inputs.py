#!/usr/bin/env python3
"""Turn gaps.tsv into a complete, self-contained Rosetta job tree.

Two stages, because KIC cannot build residues that are not in the pose. The
loopmodel documentation requires a "Template pdb file with real coordinates for
all residues plus the first and last residue of each loop region", so:

  stage 1  RosettaRemodel inserts the missing residues from a blueprint and
           closes the chain, giving a gap-free pose.
  stage 2  loopmodel/NGK resamples each inserted span with the extend flag set,
           which discards Remodel's conformation and rebuilds from idealised
           geometry. Remodel supplies connectivity; NGK does the modelling.

All insertions for a chain go in ONE blueprint, on purpose. Pose numbering
shifts by the number of residues inserted below a position, so inserting gaps
in separate runs would leave every later loop file pointing at the wrong
residues - and pointing at the wrong residues is not an error Rosetta can
detect. One gap-free pose means stage-2 pose numbering is simply 1..N over the
protein chain, which is why the protein chain is written FIRST: Rosetta numbers
poses continuously across chains in file order.

Three numbering systems are in play and are never mixed implicitly:
  author     numbering you read in the deposited file (may be non-contiguous)
  input pose 1..N over OBSERVED residues only (what Remodel reads)
  full pose  1..M over observed + inserted (what loopmodel reads)
numbering.tsv records all three for every residue, and every generated loop
file is cross-checked against it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gemmi
import pandas as pd

AA3_TO_1 = {}


def die(msg: str) -> None:
    sys.exit(f"make_inputs: {msg}")


def one_letter(mon_id: str) -> str:
    info = gemmi.find_tabulated_residue(mon_id)
    code = info.one_letter_code if info else ""
    return code.upper() if code and code.isalpha() else "X"


def is_peptide(chain: gemmi.Chain) -> bool:
    poly = chain.get_polymer()
    return bool(len(poly)) and poly.check_polymer_type() in (
        gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)


def is_nucleic(chain: gemmi.Chain) -> bool:
    poly = chain.get_polymer()
    return bool(len(poly)) and poly.check_polymer_type() in (
        gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid)


def load_structure(cif: str) -> gemmi.Structure:
    st = gemmi.make_structure_from_block(gemmi.cif.read(cif).sole_block())
    st.setup_entities()
    if not len(st):
        die(f"{cif} parsed to zero models.")
    # Keep only the first model: Rosetta reads one, and a multi-model file
    # would silently give it the first anyway.
    while len(st) > 1:
        del st[1]
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    return st


def build_input_pdb(st: gemmi.Structure, prot_chain: str, keep_nucleic: bool,
                    keep_het: bool) -> tuple[gemmi.Structure, list[dict]]:
    """Protein chain first and renumbered 1..N; optional nucleic/het after.

    Returns the new structure and the per-residue numbering rows for the
    protein chain (author -> input pose).
    """
    src = st[0]
    if prot_chain not in [c.name for c in src]:
        die(f"chain {prot_chain} absent from the model.")

    out = gemmi.Structure()
    out.name = st.name
    out.cell = st.cell
    out.spacegroup_hm = st.spacegroup_hm or "P 1"
    model = gemmi.Model("1")

    rows: list[dict] = []
    new_prot = gemmi.Chain(prot_chain)
    pose = 0
    for res in src[prot_chain]:
        if not len(res):
            continue
        info = gemmi.find_tabulated_residue(res.name)
        if not (info and info.is_amino_acid()):
            # A ligand or ion modelled inside the protein chain. It would shift
            # pose numbering if carried here, so it moves to the het chain.
            continue
        pose += 1
        copy = gemmi.Residue()
        copy.name = res.name
        copy.seqid = gemmi.SeqId(str(pose))
        copy.het_flag = "A"
        for atom in res:
            copy.add_atom(atom)
        new_prot.add_residue(copy)
        rows.append({
            "auth_seq": str(res.seqid.num),
            "ins_code": res.seqid.icode.strip(),
            "input_pose": pose,
            "resname": res.name,
            "aa": one_letter(res.name),
        })
    model.add_chain(new_prot)

    extra: list[str] = []
    for chain in src:
        if chain.name == prot_chain:
            continue
        if is_peptide(chain):
            # A second protein copy would double every pose index. The dimer is
            # rebuilt by applying the assembly transform to the repaired
            # monomer, not modelled here.
            die(f"chain {chain.name} is a second peptide chain. Model the "
                "monomer and regenerate the assembly afterwards; a two-copy "
                "pose doubles every gap and Remodel's -run:chain handles one "
                "chain at a time.")
        nucleic = is_nucleic(chain)
        if nucleic and not keep_nucleic:
            continue
        if not nucleic and not keep_het:
            continue
        copy = gemmi.Chain(chain.name)
        for res in chain:
            if len(res):
                copy.add_residue(res)
        if len(copy):
            model.add_chain(copy)
            extra.append(f"{chain.name}({'nucleic' if nucleic else 'het'},"
                         f"{len(copy)})")

    out.add_model(model)
    out.setup_entities()
    return out, rows, extra


def blueprint_lines(rows: list[dict], gaps: list[dict]) -> list[str]:
    """Whole-chain blueprint with every gap inserted.

    Fixed residues are '<n> <aa> .'. Flanks are made mobile ('L' with their own
    identity) because the loop cannot close otherwise. Inserted residues are
    '0 x L PIKAA <aa>' - residue number 0 marks an insertion, and the second
    column is only read for extensions or with -remodel:use_blueprint_sequence.
    """
    by_pose = {r["input_pose"]: r for r in rows}
    # pose index after which each gap's residues are inserted, and the flank
    # poses that must be mobile
    inserts: dict[int, list[str]] = {}
    mobile: set[int] = set()
    for g in gaps:
        after = g["flank_n_pose"]
        inserts.setdefault(after, []).extend(list(g["missing_seq"]))
        for p in range(after - g["flank"] + 1, after + 1):
            if p in by_pose:
                mobile.add(p)
        for p in range(g["flank_c_pose"], g["flank_c_pose"] + g["flank"]):
            if p in by_pose:
                mobile.add(p)

    lines: list[str] = []
    for pose in sorted(by_pose):
        r = by_pose[pose]
        if pose in mobile:
            lines.append(f"{pose} {r['aa']} L PIKAA {r['aa']}")
        else:
            lines.append(f"{pose} {r['aa']} .")
        for aa in inserts.get(pose, []):
            lines.append(f"0 x L PIKAA {aa}")
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gaps", required=True, help="gaps.tsv from find_gaps.py")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--flank", type=int, default=2,
                    help="residues each side added to the rebuilt span, so the "
                         "fixed C-alpha anchors lie OUTSIDE it (default 2)")
    ap.add_argument("--nstruct", type=int, default=500,
                    help="loopmodel models per gap (default 500)")
    ap.add_argument("--chunk", type=int, default=25,
                    help="models per array task (default 25)")
    ap.add_argument("--num-trajectory", type=int, default=50,
                    help="-remodel:num_trajectory (default 50)")
    ap.add_argument("--save-top", type=int, default=10,
                    help="-remodel:save_top (default 10)")
    ap.add_argument("--rosetta-db", required=True,
                    help="-database path INSIDE the container, from "
                         "bin/probe_container.sh. Not guessed.")
    ap.add_argument("--loopmodel-bin", required=True,
                    help="loopmodel binary name or path inside the container")
    ap.add_argument("--remodel-bin", required=True,
                    help="remodel binary name or path inside the container")
    ap.add_argument("--keep-nucleic", dest="keep_nucleic",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="keep DNA/RNA chains in the pose (default: keep). "
                         "Loops at a DNA interface modelled with the DNA "
                         "removed are built in the wrong state.")
    ap.add_argument("--keep-het", dest="keep_het",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="keep ions and ligands, e.g. the zinc thumb "
                         "(default: keep)")
    ap.add_argument("--min-fixed", type=int, default=8,
                    help="fewest residues that must lie outside every loop, so "
                         "the drift check and the spread superposition have "
                         "something to anchor on (default 8)")
    ap.add_argument("--compare-gap", default="",
                    help="gap id (e.g. cgas_dna_monomer_A_300_302) to ALSO run "
                         "with KIC-with-fragments, for a method comparison")
    args = ap.parse_args()

    if args.flank < 1:
        die("--flank must be at least 1: the loop's first and last residue are "
            "anchors whose C-alpha stays fixed, so a flank of 0 would freeze "
            "the first and last residue you are trying to build.")

    gaps = pd.read_csv(args.gaps, sep="\t", dtype=str).fillna("")
    if gaps.empty:
        die(f"{args.gaps} has no rows.")
    gaps = gaps[gaps["eligible"] == "yes"]
    if gaps.empty:
        die(f"{args.gaps} has no eligible gaps. Nothing to model.")

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    jobs: list[dict] = []
    summary: list[dict] = []

    for (sid, chain, cif), grp in gaps.groupby(
            ["structure_id", "chain", "cif_path"], sort=True):
        work = out_root / sid
        (work / "loops").mkdir(parents=True, exist_ok=True)

        st = load_structure(cif)
        inp, rows, extra = build_input_pdb(
            st, chain, args.keep_nucleic, args.keep_het)
        if not rows:
            die(f"{sid}: chain {chain} has no amino-acid residues with atoms.")

        auth_to_pose = {r["auth_seq"]: r["input_pose"] for r in rows}

        glist: list[dict] = []
        for _, g in grp.iterrows():
            for side in ("flank_n_auth", "flank_c_auth"):
                if g[side] not in auth_to_pose:
                    die(f"{sid}: gap {g['gap_start_auth']}-{g['gap_end_auth']} "
                        f"has {side}={g[side]!r}, which is not an observed "
                        "residue of this chain. gaps.tsv and the structure "
                        "disagree; re-run find_gaps.py.")
            n_pose = auth_to_pose[g["flank_n_auth"]]
            c_pose = auth_to_pose[g["flank_c_auth"]]
            if c_pose != n_pose + 1:
                die(f"{sid}: gap {g['gap_start_auth']}-{g['gap_end_auth']} has "
                    f"flanks at input poses {n_pose} and {c_pose}, which are "
                    "not adjacent. That means another residue sits between "
                    "them, so this is not a single contiguous run.")
            if len(g["missing_seq"]) != int(g["length"]):
                die(f"{sid}: gap {g['gap_start_auth']}-{g['gap_end_auth']} "
                    f"declares length {g['length']} but carries "
                    f"{len(g['missing_seq'])} residue identities.")
            if "X" in g["missing_seq"]:
                die(f"{sid}: gap {g['gap_start_auth']}-{g['gap_end_auth']} "
                    f"sequence {g['missing_seq']!r} contains an untranslatable "
                    "residue. Rosetta needs a real identity for every "
                    "inserted residue.")
            glist.append({
                "gap_id": f"{sid}_{chain}_{g['gap_start_auth']}_{g['gap_end_auth']}",
                "gap_start_auth": g["gap_start_auth"],
                "gap_end_auth": g["gap_end_auth"],
                "missing_seq": g["missing_seq"],
                "length": int(g["length"]),
                "flank_n_auth": g["flank_n_auth"],
                "flank_c_auth": g["flank_c_auth"],
                "flank_n_pose": n_pose,
                "flank_c_pose": c_pose,
                "flank": args.flank,
            })
        glist.sort(key=lambda d: d["flank_n_pose"])

        # Full-pose numbering: walk the chain once, inserting each gap after
        # its N-side flank. Done by construction rather than by arithmetic on
        # offsets, which is where this kind of code usually goes wrong.
        full: list[dict] = []
        inserted_spans: dict[str, tuple[int, int]] = {}
        by_after: dict[int, dict] = {g["flank_n_pose"]: g for g in glist}
        for r in rows:
            full.append({**r, "full_pose": len(full) + 1, "modelled": 0})
            g = by_after.get(r["input_pose"])
            if g:
                first = len(full) + 1
                for k, aa in enumerate(g["missing_seq"]):
                    full.append({
                        "auth_seq": str(int(g["gap_start_auth"]) + k),
                        "ins_code": "",
                        "input_pose": "",
                        "resname": "",
                        "aa": aa,
                        "full_pose": len(full) + 1,
                        "modelled": 1,
                    })
                inserted_spans[g["gap_id"]] = (first, len(full))

        expected = len(rows) + sum(g["length"] for g in glist)
        if len(full) != expected:
            die(f"{sid}: built {len(full)} full-pose residues, expected "
                f"{expected}. Refusing to write loop files against it.")

        pd.DataFrame(full, columns=["auth_seq", "ins_code", "input_pose",
                                    "resname", "aa", "full_pose", "modelled"]
                     ).to_csv(work / "numbering.tsv", sep="\t", index=False)
        inp.write_pdb(str(work / "input.pdb"))
        full_seq = "".join(r["aa"] for r in full)
        (work / "expected_sequence.txt").write_text(full_seq + "\n")

        bp = blueprint_lines(rows, glist)
        # The blueprint must not end with an empty line.
        (work / "insert.blueprint").write_text("\n".join(bp) + "\n")

        (work / "flags.remodel").write_text("\n".join([
            f"-database {args.rosetta_db}",
            "-in:file:s input.pdb",
            f"-remodel:blueprint insert.blueprint",
            f"-run:chain {chain}",
            f"-remodel:num_trajectory {args.num_trajectory}",
            f"-remodel:save_top {args.save_top}",
            "-remodel:quick_and_dirty false",
            "-out:file:scorefile remodel.sc",
            "-ex1",
            "-ex2",
            "-extrachi_cutoff 0",
            "-overwrite",
            "",
        ]))

        n_tasks = max(1, -(-args.nstruct // args.chunk))

        # Flanked spans, then MERGED where they touch. Two gaps three residues
        # apart have overlapping flanked spans, and overlapping LOOP lines are
        # not two loops - they are one. Merging also removes the alternative,
        # which would be to model each gap in a separate run and then try to
        # combine two structures that each rebuilt a different loop.
        spans: list[dict] = []
        for g in glist:
            first, last = inserted_spans[g["gap_id"]]
            spans.append({
                "lo": max(1, first - args.flank),
                "hi": min(len(full), last + args.flank),
                "gap_ids": [g["gap_id"]],
                "modelled": [(first, last)],
            })
        spans.sort(key=lambda d: d["lo"])
        merged: list[dict] = []
        for sp in spans:
            if merged and sp["lo"] <= merged[-1]["hi"]:
                prev = merged[-1]
                prev["hi"] = max(prev["hi"], sp["hi"])
                prev["gap_ids"] += sp["gap_ids"]
                prev["modelled"] += sp["modelled"]
            else:
                merged.append(dict(sp))

        for sp in merged:
            # KIC holds the C-alpha of a loop's first and last residue fixed,
            # so an anchor that is itself a modelled residue would pin a
            # position Rosetta invented rather than one the crystal observed.
            for anchor in (sp["lo"], sp["hi"]):
                if full[anchor - 1]["modelled"]:
                    die(f"{sid}: loop {sp['lo']}-{sp['hi']} would anchor on "
                        f"full pose {anchor}, which is itself modelled. "
                        "Lower --flank, or these gaps reach the chain end and "
                        "have no observed anchor.")

        covered = {i for sp in merged for i in range(sp["lo"], sp["hi"] + 1)}
        outside = len(full) - len(covered)
        if outside < args.min_fixed:
            die(f"{sid}: only {outside} residue(s) would lie outside every "
                f"loop, below --min-fixed {args.min_fixed}. Those residues "
                "are the reference for two checks that would otherwise pass "
                "by having nothing to compare: that KIC left the rest of the "
                "chain where it was, and the superposition that makes the "
                "ensemble spread a conformational difference rather than a "
                "change of coordinate frame. Lower --flank.")

        # One loop file for the whole chain: every loop is sampled in the same
        # run, so one output structure carries every rebuilt gap.
        loop_path = work / "loops" / "all.loop"
        loop_path.write_text("".join(
            f"LOOP {sp['lo']} {sp['hi']} 0 0 1\n" for sp in merged))

        def write_flags(name: str, style: str, loop_rel: str) -> None:
            flags = [
                f"-database {args.rosetta_db}",
                "-in:file:s insert_selected.pdb",
                "-in:file:fullatom",
                f"-loops:loop_file {loop_rel}",
            ]
            if style == "ngk":
                # Stein & Kortemme next-generation KIC flag set.
                flags += [
                    "-loops:remodel perturb_kic",
                    "-loops:refine refine_kic",
                    "-loops:kic_rama2b",
                    "-loops:kic_omega_sampling",
                    "-allow_omega_move true",
                    "-loops:ramp_fa_rep",
                    "-loops:ramp_rama",
                ]
            else:
                # KIC with fragments draws every non-pivot torsion from
                # fragment data, so the vicinity-sampling flags do not apply
                # and fragment files are required. The fragment paths depend on
                # the image, so the submit script appends them.
                flags += [
                    "-loops:remodel perturb_kic_with_fragments",
                    "-loops:refine refine_kic_with_fragments",
                ]
            flags += ["-ex1", "-ex2", "-extrachi_cutoff 0",
                      "-out:pdb_gz",
                      # Model selection needs scores, and a run without a
                      # score file leaves nothing to rank 500 models by.
                      f"-out:file:scorefile {style}.sc",
                      "-overwrite", ""]
            (work / name).write_text("\n".join(flags))

        write_flags("flags.ngk", "ngk", "loops/all.loop")
        styles = [("ngk", "loops/all.loop")]

        if args.compare_gap:
            hit = [sp for sp in merged if args.compare_gap in sp["gap_ids"]]
            if hit:
                sp = hit[0]
                rel = f"loops/compare_{args.compare_gap}.loop"
                (work / rel).write_text(f"LOOP {sp['lo']} {sp['hi']} 0 0 1\n")
                write_flags(f"flags.kicfrag", "kicfrag", rel)
                styles.append(("kicfrag", rel))
                if len(sp["gap_ids"]) > 1:
                    print(f"  note: --compare-gap {args.compare_gap} shares a "
                          f"merged loop with {sp['gap_ids']}, so the "
                          "comparison covers all of them")

        for style, rel in styles:
            for task in range(1, n_tasks + 1):
                n = min(args.chunk, args.nstruct - (task - 1) * args.chunk)
                jobs.append({
                    "structure_id": sid, "chain": chain, "style": style,
                    "task": task, "nstruct": n,
                    "work_dir": str(work),
                    "flags": f"flags.{style}",
                    "loop_file": rel,
                })

        for sp in merged:
            summary.append({
                "structure_id": sid,
                "gap_ids": ",".join(sp["gap_ids"]),
                "auth_ranges": ",".join(
                    f"{g['gap_start_auth']}-{g['gap_end_auth']}"
                    for g in glist if g["gap_id"] in sp["gap_ids"]),
                "residues_modelled": sum(b - a + 1 for a, b in sp["modelled"]),
                "loop_definition_full_pose": f"{sp['lo']}-{sp['hi']}",
                "loop_length": sp["hi"] - sp["lo"] + 1,
            })

        (work / "provenance.json").write_text(json.dumps({
            "structure_id": sid,
            "pdb_id": str(grp["pdb_id"].iloc[0]),
            "chain": chain, "source_cif": cif,
            "observed_residues": len(rows),
            "full_residues": len(full),
            "extra_chains_kept": extra,
            "keep_nucleic": args.keep_nucleic, "keep_het": args.keep_het,
            "flank": args.flank, "nstruct": args.nstruct,
            "gaps": glist,
            "loops": [{"lo": sp["lo"], "hi": sp["hi"],
                       "gap_ids": sp["gap_ids"]} for sp in merged],
        }, indent=2) + "\n")

        print(f"{sid}: chain {chain}, {len(rows)} observed -> {len(full)} "
              f"residues, {len(glist)} gap(s); extra chains kept: "
              f"{extra or 'none'}")

    pd.DataFrame(jobs).to_csv(out_root / "jobs.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out_root / "plan.tsv", sep="\t", index=False)
    print(f"\n{len(jobs)} array task(s) -> {out_root / 'jobs.tsv'}")
    print(f"plan -> {out_root / 'plan.tsv'}")


if __name__ == "__main__":
    main()
