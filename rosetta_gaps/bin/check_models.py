#!/usr/bin/env python3
"""Validate Rosetta output, and measure how determined each rebuilt loop is.

Two modes, one definition of "valid" shared between them:

  --stage insert   Check RosettaRemodel's candidates (1.pdb, 2.pdb, ...) and
                   write insert_selected.pdb. Candidates are taken in order and
                   the first VALID one wins, not the lowest-scoring one: stage 2
                   sets the loop file's extend flag, which discards this
                   conformation and rebuilds from idealised geometry, so all
                   this structure has to do is be complete and closed. Ranking
                   it by energy would imply a choice that has no consequence.

  --stage ngk      Check every loopmodel model, rank by total_score, and report
                   the ENSEMBLE SPREAD per gap: pairwise C-alpha RMSD over the
                   rebuilt residues across the top models, after superposing on
                   the residues that were never rebuilt.

Why the spread matters more than the score. In AICG2+ the reference structure
defines the native contacts and the equilibrium bond, angle and dihedral
values, so the loop conformation carried into the simulation BECOMES the
definition of the native state. Nothing downstream can penalise a wrong loop.
A gap whose top models agree is forced by the surrounding structure; a gap
whose top models disagree by several Angstrom is a modelling assumption, and
belongs in the limitations of any claim that rests on it.

Checks applied to every model:
  1. sequence identical to expected_sequence.txt (no residue lost, inserted
     twice, or mutated by the rebuild)
  2. no chain break: every consecutive C(i)-N(i+1) distance inside
     --bond-range, measured across the whole chain and not only at the former
     gap edges
  3. residues outside every loop definition have not moved: C-alpha RMSD to
     the stage input must stay under --fixed-tol. KIC should hold them fixed,
     so a nonzero value means the loop file or the numbering was wrong.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import shutil
import sys
from pathlib import Path

import gemmi
import pandas as pd

# A peptide bond is about 1.33 A. The window is wide enough not to fire on an
# unidealised but connected backbone, and narrow enough that an unclosed break
# (which lands in the many-Angstrom range) cannot pass.
DEFAULT_BOND_RANGE = (1.10, 1.70)


def die(msg: str) -> None:
    sys.exit(f"check_models: {msg}")


def read_model(path: Path) -> gemmi.Structure:
    if path.suffix == ".gz":
        text = gzip.decompress(path.read_bytes()).decode()
        st = gemmi.read_pdb_string(text)
    else:
        st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def protein_chain(st: gemmi.Structure, name: str) -> gemmi.Chain:
    for chain in st[0]:
        if chain.name == name:
            return chain
    die(f"chain {name} absent from {st.name}; chains present: "
        f"{[c.name for c in st[0]]}")


def ca_positions(chain: gemmi.Chain) -> list[gemmi.Position | None]:
    out = []
    for res in chain:
        atom = res.find_atom("CA", "*")
        out.append(atom.pos if atom else None)
    return out


def sequence(chain: gemmi.Chain) -> str:
    letters = []
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        code = info.one_letter_code if info else ""
        letters.append(code.upper() if code and code.isalpha() else "X")
    return "".join(letters)


def chain_breaks(chain: gemmi.Chain, lo: float, hi: float) -> list[str]:
    bad = []
    residues = list(chain)
    for i in range(len(residues) - 1):
        c = residues[i].find_atom("C", "*")
        n = residues[i + 1].find_atom("N", "*")
        if c is None or n is None:
            bad.append(f"{residues[i].seqid.num}->{residues[i+1].seqid.num}: "
                       f"missing backbone {'C' if c is None else 'N'}")
            continue
        d = c.pos.dist(n.pos)
        if not (lo <= d <= hi):
            bad.append(f"{residues[i].seqid.num}->{residues[i+1].seqid.num}: "
                       f"C-N {d:.2f} A")
    return bad


def direct_rmsd(ref: list[gemmi.Position], mob: list[gemmi.Position]) -> float:
    return math.sqrt(sum(r.dist(m) ** 2 for r, m in zip(ref, mob)) / len(ref))


def fit(ref: list[gemmi.Position], mob: list[gemmi.Position],
        eps: float = 1e-3) -> tuple[gemmi.Transform, float, str]:
    """Superpose mob onto ref, tolerating the case where they coincide.

    KIC holds every residue outside the loop fixed, so across two models of the
    same run those coordinates are usually bit-identical - and Kabsch is
    undefined for identical point sets: gemmi prints "More than 50 iterations
    needed!" and returns NaN. NaN then passes every `> tolerance` test, because
    any comparison with NaN is false, so the drift check would report no drift
    precisely when the calculation failed. Return the identity instead, and
    report a NaN from a genuinely non-identical pair rather than hiding it.
    """
    direct = direct_rmsd(ref, mob)
    if direct < eps:
        return gemmi.Transform(), direct, ""
    sup = gemmi.superpose_positions(ref, mob)
    if math.isnan(sup.rmsd):
        return gemmi.Transform(), direct, (
            f"superposition returned NaN over {len(ref)} points whose direct "
            f"RMSD is {direct:.3f} A; falling back to no superposition")
    return sup.transform, sup.rmsd, ""


def loop_spread(ranked, fixed_idx: list[int],
                loop_idx: list[int]) -> list[float]:
    """Pairwise C-alpha RMSD over loop_idx, after aligning on fixed_idx.

    Aligning first is what makes this a comparison of loop CONFORMATION rather
    than of a rigid-body offset between two output files.
    """
    dists: list[float] = []
    for (_, _, ca_a), (_, _, ca_b) in itertools.combinations(ranked, 2):
        tr, _, _ = fit([ca_a[i] for i in fixed_idx],
                       [ca_b[i] for i in fixed_idx])
        total = 0.0
        for i in loop_idx:
            v = tr.apply(ca_b[i])
            total += ca_a[i].dist(gemmi.Position(v.x, v.y, v.z)) ** 2
        d = math.sqrt(total / len(loop_idx))
        if math.isnan(d):
            return []
        dists.append(d)
    return dists


def load_context(work: Path) -> dict:
    prov = json.loads((work / "provenance.json").read_text())
    numbering = pd.read_csv(work / "numbering.tsv", sep="\t", dtype=str).fillna("")
    expected = (work / "expected_sequence.txt").read_text().strip()
    return {"prov": prov, "numbering": numbering, "expected": expected}


def loop_indices(prov: dict) -> list[int]:
    """0-based full-pose indices covered by any loop definition."""
    idx: set[int] = set()
    for lp in prov.get("loops", []):
        idx.update(range(lp["lo"] - 1, lp["hi"]))
    return sorted(idx)


def modelled_indices(prov: dict, numbering: pd.DataFrame) -> dict[str, list[int]]:
    """0-based full-pose indices of the inserted residues, per gap."""
    out: dict[str, list[int]] = {}
    for g in prov["gaps"]:
        lo, hi = int(g["gap_start_auth"]), int(g["gap_end_auth"])
        rows = numbering[(numbering["modelled"] == "1")
                         & (numbering["auth_seq"].astype(int) >= lo)
                         & (numbering["auth_seq"].astype(int) <= hi)]
        out[g["gap_id"]] = [int(v) - 1 for v in rows["full_pose"]]
    return out


def read_scores(path: Path) -> dict[str, float]:
    """description -> total_score from one Rosetta score file."""
    if not path.exists():
        return {}
    scores: dict[str, float] = {}
    header: list[str] = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "SCORE:" and parts[1] == "total_score":
            header = parts[1:]
        elif parts[0] == "SCORE:" and header:
            row = dict(zip(header, parts[1:]))
            desc = row.get("description")
            try:
                if desc:
                    scores[desc] = float(row["total_score"])
            except (KeyError, ValueError):
                continue
    return scores


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", required=True,
                    help="per-structure directory made by make_inputs.py")
    ap.add_argument("--stage", choices=("insert", "ngk"), required=True)
    ap.add_argument("--models", nargs="*", default=[],
                    help="model files. Default for insert: <work>/1.pdb.. ; "
                         "for ngk: <work>/models/*.pdb(.gz)")
    ap.add_argument("--top", type=int, default=10,
                    help="models used for the spread report (default 10)")
    ap.add_argument("--fixed-tol", type=float, default=0.10,
                    help="C-alpha RMSD tolerated outside every loop, Angstrom "
                         "(default 0.10)")
    ap.add_argument("--bond-range", nargs=2, type=float,
                    default=list(DEFAULT_BOND_RANGE),
                    help="acceptable consecutive C-N distance, Angstrom")
    ap.add_argument("--report", default="",
                    help="TSV to write (default <work>/<stage>_report.tsv)")
    args = ap.parse_args()

    work = Path(args.work_dir)
    if not (work / "provenance.json").exists():
        die(f"{work} has no provenance.json. Run make_inputs.py first.")
    ctx = load_context(work)
    prov, numbering, expected = ctx["prov"], ctx["numbering"], ctx["expected"]
    chain_name = prov["chain"]

    if args.models:
        models = [Path(m) for m in args.models]
    elif args.stage == "insert":
        models = sorted(work.glob("[0-9]*.pdb"),
                        key=lambda p: int(p.stem) if p.stem.isdigit() else 1 << 30)
    else:
        # Array tasks each run in their own runs/t<N>/ directory so that
        # Rosetta's output names, which restart at _0001 in every task, cannot
        # collide. Each task writes its own score file beside its models.
        models = sorted(work.glob("runs/*/*.pdb.gz")) \
            + sorted(work.glob("runs/*/*.pdb")) \
            + sorted((work / "models").glob("*.pdb.gz")) \
            + sorted((work / "models").glob("*.pdb"))
    if not models:
        die(f"no model files found for stage {args.stage} under {work}. "
            "Nothing was produced, or the submit script wrote elsewhere - "
            "this is a failure, not an empty result.")

    loops = loop_indices(prov)
    n_full = len(numbering)
    if args.stage == "insert":
        # Remodel's output has the full sequence; the only positions that can
        # legitimately differ from the input are the inserted residues and
        # their mobile flanks, i.e. the loop spans.
        ref_st = read_model(work / "input.pdb")
        ref_chain = protein_chain(ref_st, chain_name)
        ref_ca_input = ca_positions(ref_chain)
        # Map input-pose -> full-pose so the comparison is position-correct.
        pairs = [(int(r["input_pose"]) - 1, int(r["full_pose"]) - 1)
                 for _, r in numbering.iterrows() if r["input_pose"]]
        fixed_pairs = [(i, f) for i, f in pairs if f not in set(loops)]
    else:
        ref_st = read_model(work / "insert_selected.pdb")
        ref_ca_input = ca_positions(protein_chain(ref_st, chain_name))
        fixed_pairs = [(i, i) for i in range(n_full) if i not in set(loops)]

    rows: list[dict] = []
    good: list[tuple[Path, gemmi.Structure, list[gemmi.Position]]] = []
    for path in models:
        st = read_model(path)
        chain = protein_chain(st, chain_name)
        ca = ca_positions(chain)
        seq = sequence(chain)
        problems: list[str] = []
        if seq != expected:
            problems.append(
                f"sequence differs from expected ({len(seq)} vs "
                f"{len(expected)} residues)")
        else:
            breaks = chain_breaks(chain, *args.bond_range)
            if breaks:
                problems.append(f"{len(breaks)} backbone break(s): "
                                + "; ".join(breaks[:3]))
            if any(p is None for p in ca):
                problems.append("residue(s) without a C-alpha")
            elif fixed_pairs:
                ref_sel = [ref_ca_input[i] for i, _ in fixed_pairs]
                mob_sel = [ca[f] for _, f in fixed_pairs]
                if any(p is None for p in ref_sel):
                    problems.append("reference residue without a C-alpha")
                else:
                    _, drift, note = fit(ref_sel, mob_sel)
                    if note:
                        problems.append(note)
                    if math.isnan(drift):
                        problems.append(
                            "drift over the fixed residues is NaN: the check "
                            "did not run, so do not read this model as clean")
                    elif drift > args.fixed_tol:
                        problems.append(
                            f"residues outside every loop moved: C-alpha RMSD "
                            f"{drift:.2f} A > {args.fixed_tol}")
        rows.append({"model": path, "valid": "no" if problems else "yes",
                     "problems": "; ".join(problems)})
        if not problems:
            good.append((path, st, ca))

    report = Path(args.report) if args.report else work / f"{args.stage}_report.tsv"

    if args.stage == "insert":
        flat = [{**r, "model": r["model"].name} for r in rows]
        pd.DataFrame(flat).to_csv(report, sep="\t", index=False)
        if not good:
            print(pd.DataFrame(flat).to_string(index=False), file=sys.stderr)
            die("no Remodel candidate passed. Nothing is selected, so stage 2 "
                "has no input. The report above says why for each candidate.")
        chosen = good[0][0]
        dest = work / "insert_selected.pdb"
        if chosen.suffix == ".gz":
            dest.write_bytes(gzip.decompress(chosen.read_bytes()))
        else:
            shutil.copyfile(chosen, dest)
        print(f"{len(good)}/{len(models)} candidate(s) valid; selected "
              f"{chosen.name} -> {dest}")
        print(f"report -> {report}")
        return

    # ---- ngk: rank, then measure spread over the rebuilt residues ----
    # Scores come from the score file BESIDE each model, not from one global
    # file: two array tasks both produce a model described as "..._0001", and a
    # single lookup table would silently give one of them the other's score.
    cache: dict[Path, dict[str, float]] = {}

    def score_of(p: Path) -> float:
        d = p.parent
        if d not in cache:
            sc = d / "ngk.sc"
            if not sc.exists():
                candidates = sorted(d.glob("*.sc"))
                sc = candidates[0] if candidates else sc
            cache[d] = read_scores(sc)
        stem = p.name
        for suffix in (".pdb.gz", ".pdb"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
        return cache[d].get(stem, float("inf"))

    def label(p: Path) -> str:
        try:
            return str(p.relative_to(work))
        except ValueError:
            return p.name

    good.sort(key=lambda t: score_of(t[0]))
    ranked = good[: args.top]
    missing_scores = sum(1 for p, _, _ in good if score_of(p) == float("inf"))
    if missing_scores:
        print(f"warning: {missing_scores} valid model(s) have no total_score "
              "in the score file beside them and sort last. Ranking is "
              "incomplete - check that each run directory holds its own "
              "*.sc.", file=sys.stderr)

    per_model = pd.DataFrame([{
        "model": label(r["model"]),
        "valid": r["valid"],
        "problems": r["problems"],
        "total_score": ("" if score_of(r["model"]) == float("inf")
                        else f"{score_of(r['model']):.3f}"),
    } for r in rows])
    per_model.to_csv(report, sep="\t", index=False)

    gaps = modelled_indices(prov, numbering)
    fixed_idx = [f for _, f in fixed_pairs]
    spread_rows: list[dict] = []
    for gap_id, idx in gaps.items():
        if len(ranked) < 2:
            spread_rows.append({"gap_id": gap_id, "n_models": len(ranked),
                                "mean_pairwise_rmsd": "",
                                "max_pairwise_rmsd": "",
                                "verdict": "too few valid models to judge"})
            continue
        dists = loop_spread(ranked, fixed_idx, idx)
        if not dists:
            spread_rows.append({"gap_id": gap_id, "n_models": len(ranked),
                                "mean_pairwise_rmsd": "",
                                "max_pairwise_rmsd": "",
                                "verdict": "spread could not be computed"})
            continue
        mean_d, max_d = sum(dists) / len(dists), max(dists)
        verdict = ("tightly determined" if max_d < 1.0 else
                   "moderately determined" if max_d < 2.5 else
                   "NOT determined by the data")
        spread_rows.append({"gap_id": gap_id, "n_models": len(ranked),
                            "mean_pairwise_rmsd": f"{mean_d:.2f}",
                            "max_pairwise_rmsd": f"{max_d:.2f}",
                            "verdict": verdict})

    spread = work / "spread_report.tsv"
    pd.DataFrame(spread_rows).to_csv(spread, sep="\t", index=False)

    print(f"{len(good)}/{len(models)} model(s) valid")
    if ranked:
        best = ranked[0][0]
        shown = (f"{score_of(best):.3f}" if score_of(best) != float("inf")
                 else "no score")
        print(f"best by total_score: {label(best)} ({shown})")
    print(f"\nensemble spread over the top {len(ranked)} valid model(s), "
          "C-alpha RMSD of the rebuilt residues after superposing on the rest:")
    print(pd.DataFrame(spread_rows).to_string(index=False))
    print(f"\nper-model report -> {report}")
    print(f"spread report     -> {spread}")
    print("\nA gap reported as NOT determined is a modelling assumption, not a "
          "result: in AICG2+ whichever conformation you promote becomes the "
          "native state and nothing downstream can penalise it.")


if __name__ == "__main__":
    main()
