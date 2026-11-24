#!/usr/bin/env python3
"""Prepare cGAS structures for CG modelling, reverting crystallisation mutations to WT.

Why this step exists
--------------------
Several cGAS structures in the PDB carry substitutions introduced to make the
protein crystallise. They sit in or near the DNA-binding interface, which is
the region this pipeline measures, so leaving them in place confounds exactly
the comparison being made. This step reverts them and records the conversion.

How the mutations are identified
--------------------------------
NOT from a hardcoded table. The structure's own sequence is aligned to a
user-supplied reference sequence (the canonical human cGAS) and every
difference is reported with its author residue number. That is verifiable,
and it generalises to any structure passed in. If the mmCIF carries
_struct_ref_seq_dif records, they are read as an independent cross-check and
any disagreement with the alignment is reported rather than resolved silently.

Why a rename is not enough
--------------------------
AICG2+ is a structure-based model: its native contact map is derived from
all-atom coordinates, sidechains included. Renaming LYS to ASN would leave the
crystallographic lysine sidechain in place, so contacts would still be
computed from atoms that do not exist in the wild-type residue.

  --sidechain-policy strip  (default)
      Keep backbone N, CA, C, O and CB at reverted positions; drop every
      atom beyond CB. No contact can then be derived from a sidechain that
      the wild-type residue does not have. The cost is that genuine
      wild-type sidechain contacts at those positions are lost too, which is
      conservative in the right direction: a missing contact weakens the
      native-state bias, an invented one would bias the result.
      Reverted positions are reported so their effect can be bounded.

  --sidechain-policy label-only
      Rename and keep all atoms. Faster, and wrong in the way described
      above. Available for diagnosis, not recommended for production.

Full sidechain rebuilding would need a side-chain packer, which is a separate
tool and a separate decision; this script does not pretend to do it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BACKBONE_PLUS_CB = {"N", "CA", "C", "O", "CB", "OXT"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if k not in ("MSE", "SEC", "PYL")}


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def read_reference(path: str) -> str:
    seq, seen = [], 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                seen += 1
                if seen > 1:
                    break
                continue
            seq.append(line)
    if not seq:
        die(f"{path} contains no sequence")
    s = "".join(seq).upper()
    bad = set(s) - set(ONE_TO_THREE)
    if bad:
        die(f"reference sequence has non-standard residue letter(s): {sorted(bad)}")
    return s


def chain_residues(st, chain_name: str):
    """Observed polymer residues of one chain, in order, with author numbering."""
    model = st[0]
    ch = None
    for c in model:
        if c.name == chain_name:
            ch = c
            break
    if ch is None:
        have = ",".join(c.name for c in model)
        die(f"chain {chain_name!r} not in model; chains present: {have}")
    out = []
    for res in ch:
        if res.name not in THREE_TO_ONE:
            continue
        out.append(res)
    if not out:
        die(f"chain {chain_name!r} has no standard amino-acid residues")
    return ch, out


def align(ref: str, obs: str):
    from Bio import Align
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    # end_gap_score is the stable umbrella attribute; the per-side names were
    # renamed in Biopython 1.88 and setting them by the old names warns.
    al.end_gap_score = 0.0
    try:
        best = al.align(ref, obs)[0]
    except IndexError:
        die("alignment produced no result; is the reference the right protein?")
    return best


def cif_seq_dif(doc) -> list[dict]:
    """Read _struct_ref_seq_dif as an independent cross-check, if present."""
    out = []
    try:
        block = doc.sole_block()
        tab = block.find("_struct_ref_seq_dif.",
                         ["pdbx_pdb_strand_id", "pdbx_auth_seq_num",
                          "db_mon_id", "mon_id", "details"])
    except Exception:
        return out
    for row in tab:
        try:
            out.append({
                "chain": str(row[0]), "auth_seq": str(row[1]),
                "db_residue": str(row[2]), "pdb_residue": str(row[3]),
                "details": str(row[4]),
            })
        except Exception:
            continue
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True,
                    help="TSV: structure_id, pdb_id, chain, cif_path[, revert]")
    ap.add_argument("--reference", required=True,
                    help="FASTA of the canonical reference protein sequence")
    ap.add_argument("--sidechain-policy", choices=["strip", "label-only"],
                    default="strip")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-substitutions", type=int, default=10,
                    help="more differences than this is treated as the wrong "
                         "reference or the wrong chain, and fails")
    args = ap.parse_args()

    import gemmi
    import pandas as pd

    ref = read_reference(args.reference)
    out = Path(args.out_dir)
    (out / "structures").mkdir(parents=True, exist_ok=True)

    ss = pd.read_csv(args.samplesheet, sep="\t", dtype=str).fillna("")
    for col in ("structure_id", "chain", "cif_path"):
        if col not in ss.columns:
            die(f"samplesheet is missing required column {col!r}")
    if ss["structure_id"].duplicated().any():
        dup = ss.loc[ss["structure_id"].duplicated(), "structure_id"].tolist()
        die(f"duplicate structure_id(s): {dup}")

    records = []
    for _, row in ss.iterrows():
        sid, chain, cif = row["structure_id"], row["chain"], row["cif_path"]
        revert = (row.get("revert") or "auto").strip() or "auto"
        if not cif:
            die(f"{sid}: cif_path is empty. This script does not download "
                "structures; fetch them in the pipeline's fetch step (which "
                "needs network access) and pass the local path here.")
        if not Path(cif).exists():
            die(f"{sid}: cif_path does not exist: {cif}")

        doc = gemmi.cif.read(cif)
        st = gemmi.make_structure_from_block(doc.sole_block())
        st.setup_entities()
        ch, residues = chain_residues(st, chain)
        obs = "".join(THREE_TO_ONE[r.name] for r in residues)

        aln = align(ref, obs)
        ti, qi = aln.aligned  # paired blocks of (start, end) in ref / obs
        subs, ins, dele = [], [], []
        paired_ref, paired_obs = set(), set()
        for (rs, re_), (qs, qe) in zip(ti, qi):
            for k in range(re_ - rs):
                ri, oi = rs + k, qs + k
                paired_ref.add(ri)
                paired_obs.add(oi)
                if ref[ri] != obs[oi]:
                    res = residues[oi]
                    # ints here come from numpy alignment arrays; cast so the
                    # record is JSON-serialisable.
                    subs.append({
                        "auth_seq": int(res.seqid.num),
                        "icode": res.seqid.icode.strip(),
                        "ref_residue": ref[ri],
                        "pdb_residue": obs[oi],
                        "ref_index_1based": int(ri) + 1,
                    })
        ins = [{"auth_seq": int(residues[i].seqid.num), "pdb_residue": obs[i]}
               for i in range(len(obs)) if i not in paired_obs]
        dele_idx = [i for i in range(len(ref)) if i not in paired_ref]
        dele = {"n_reference_residues_absent": len(dele_idx)}

        if len(subs) > args.max_substitutions:
            die(f"{sid}: alignment to the reference shows {len(subs)} "
                f"substitutions, above --max-substitutions "
                f"{args.max_substitutions}. That is far more than "
                "crystallisation engineering would explain - check that "
                f"chain {chain!r} is the right chain and that the reference "
                "is the right protein and species.")

        xcheck = [d for d in cif_seq_dif(doc) if d["chain"] in ("", chain)]

        # Apply the reversion
        if revert.lower() == "none":
            to_revert = []
        elif revert.lower() == "auto":
            to_revert = list(subs)
        else:
            wanted = {w.strip().upper() for w in revert.split(";") if w.strip()}
            index = {f"{s['pdb_residue']}{s['auth_seq']}{s['ref_residue']}": s
                     for s in subs}
            unknown = wanted - set(index)
            if unknown:
                die(f"{sid}: revert spec names {sorted(unknown)}, which the "
                    f"alignment does not show as substitutions. Observed: "
                    f"{sorted(index)}")
            to_revert = [index[w] for w in sorted(wanted)]

        applied, atoms_dropped = [], 0
        for s in to_revert:
            target = None
            for res in residues:
                if res.seqid.num == s["auth_seq"]:
                    target = res
                    break
            if target is None:
                die(f"{sid}: cannot locate author residue {s['auth_seq']}")
            old_name, old_atoms = target.name, len(target)
            target.name = ONE_TO_THREE[s["ref_residue"]]
            if args.sidechain_policy == "strip":
                keep = [a for a in target if a.name in BACKBONE_PLUS_CB]
                if s["ref_residue"] == "G":
                    keep = [a for a in keep if a.name != "CB"]
                removed = old_atoms - len(keep)
                kept_copy = [a.clone() for a in keep]
                for _ in range(len(target)):
                    del target[0]
                for a in kept_copy:
                    target.add_atom(a)
                atoms_dropped += removed
            applied.append({
                "auth_seq": int(s["auth_seq"]),
                "from": f"{old_name}({s['pdb_residue']})",
                "to": f"{target.name}({s['ref_residue']})",
                "atoms_after": int(len(target)),
            })

        st.setup_entities()
        cif_out = out / "structures" / f"{sid}.wt.cif"
        st.make_mmcif_document().write_file(str(cif_out))

        rec = {
            "structure_id": sid,
            "pdb_id": row.get("pdb_id", ""),
            "chain": chain,
            "source_cif": str(cif),
            "prepared_cif": str(cif_out),
            "reference_fasta": args.reference,
            "reference_length": len(ref),
            "observed_residues": int(len(obs)),
            "sidechain_policy": args.sidechain_policy,
            "revert_spec": revert,
            "differences_from_reference": {
                "substitutions": subs, "insertions": ins, "deletions": dele,
            },
            "reverted": applied,
            "atoms_dropped": int(atoms_dropped),
            "struct_ref_seq_dif_crosscheck": xcheck,
            "crosscheck_agrees": None,
        }
        if xcheck:
            aln_set = {(str(s["auth_seq"]), s["pdb_residue"]) for s in subs}
            cif_set = {(d["auth_seq"], THREE_TO_ONE.get(d["pdb_residue"], "?"))
                       for d in xcheck
                       if "MUTATION" in d["details"].upper()}
            rec["crosscheck_agrees"] = bool(cif_set <= aln_set) if cif_set else None
            if cif_set and not cif_set <= aln_set:
                sys.stderr.write(
                    f"WARNING {sid}: _struct_ref_seq_dif lists engineered "
                    f"mutation(s) {sorted(cif_set - aln_set)} that the alignment "
                    "to your reference does not show as differences. Check the "
                    "reference isoform.\n")
        records.append(rec)

        print(f"[prepare_structures] {sid} chain={chain} "
              f"obs={len(obs)}/{len(ref)} subs={len(subs)} "
              f"reverted={len(applied)} atoms_dropped={atoms_dropped} "
              f"policy={args.sidechain_policy}"
              + ("" if applied else "  (no crystallisation mutations found)"))

    with open(out / "structure_preparation.json", "w") as fh:
        json.dump({"sidechain_policy": args.sidechain_policy,
                   "reference_fasta": args.reference,
                   "structures": records}, fh, indent=2)

    rows = []
    for r in records:
        for a in r["reverted"]:
            rows.append({"structure_id": r["structure_id"],
                         "pdb_id": r["pdb_id"], "chain": r["chain"],
                         "auth_seq": a["auth_seq"], "from": a["from"],
                         "to": a["to"], "policy": r["sidechain_policy"]})
    import pandas as pd  # noqa: F811
    pd.DataFrame(rows, columns=["structure_id", "pdb_id", "chain", "auth_seq",
                                "from", "to", "policy"]
                 ).to_csv(out / "wt_reversions.tsv", sep="\t", index=False)
    print(f"[prepare_structures] wrote {out/'structure_preparation.json'} "
          f"and {out/'wt_reversions.tsv'} ({len(rows)} reversion(s))")


if __name__ == "__main__":
    main()
