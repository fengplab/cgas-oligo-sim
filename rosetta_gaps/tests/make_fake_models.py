#!/usr/bin/env python3
"""Fabricate Rosetta-shaped output for check_models.py, from the real input.

Observed residues are copied verbatim out of <work>/input.pdb, so the residues
that KIC holds fixed are bit-identical to the stage input - which is what the
drift check compares against, and also the case that makes Kabsch superposition
degenerate. Only the inserted residues are synthesised: their C-alphas sit on
an arc between the flanking C-alphas, with the bulge height solved so the
spacing comes out near 3.80 A, and N/C/O follow the local chain direction so
the peptide bonds land in range by construction rather than by luck.

The backbone is not physical. These fixtures exercise the CHECKS: a model with
a planted break, a lost residue or a shifted fixed region must be rejected, and
a set of models differing only inside the loop must yield a spread.

Modes:
  insert   candidates 1.pdb..n.pdb, optionally carrying one planted defect
  ngk      models/ngk_NNNN.pdb plus ngk.sc, loop C-alphas displaced by a
           tapered random offset of the given amplitude
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

import gemmi
import pandas as pd

THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
         "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
         "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
         "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
STEP = 3.80


def vsub(a, b): return [a[i] - b[i] for i in range(3)]
def vadd(a, b): return [a[i] + b[i] for i in range(3)]
def vmul(a, s): return [c * s for c in a]
def vlen(a): return math.sqrt(sum(c * c for c in a))
def vnorm(a):
    n = vlen(a) or 1.0
    return [c / n for c in a]


def perpendicular(u):
    """Any unit vector orthogonal to u."""
    seed = [0.0, 0.0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0.0, 0.0]
    v = vsub(seed, vmul(u, sum(seed[i] * u[i] for i in range(3))))
    return vnorm(v)


def arc(p0, p1, k, h, normal):
    """k interior points from p0 to p1, bulged by h along normal."""
    pts = []
    for j in range(1, k + 1):
        f = j / (k + 1)
        base = vadd(p0, vmul(vsub(p1, p0), f))
        pts.append(vadd(base, vmul(normal, h * math.sin(math.pi * f))))
    return pts


def solve_bulge(p0, p1, k, normal):
    """Bulge height whose polyline spacing is closest to STEP.

    The chord across a gap is shorter than (k+1) * STEP, so a straight
    interpolation would compress every bond. Bisection on h is enough: path
    length grows monotonically with h.
    """
    def worst(h):
        pts = [p0] + arc(p0, p1, k, h, normal) + [p1]
        return max(abs(vlen(vsub(pts[i + 1], pts[i])) - STEP)
                   for i in range(len(pts) - 1))
    # A scan, not a bisection: worst() is not monotonic in h, so a gradient
    # step can settle on the wrong side and leave a bond out of range.
    best_h, best = 0.0, worst(0.0)
    h = 0.0
    while h <= 20.0:
        w = worst(h)
        if w < best:
            best_h, best = h, w
        h += 0.01
    return best_h


def residue(name, num, ca, direction):
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(str(num))
    res.het_flag = "A"
    u = vnorm(direction)
    for atom_name, elem, offset in (("N", "N", -0.95), ("CA", "C", 0.0),
                                    ("C", "C", 1.52), ("O", "O", 2.10)):
        at = gemmi.Atom()
        at.name, at.element = atom_name, gemmi.Element(elem)
        at.occ, at.b_iso = 1.0, 20.0
        p = vadd(ca, vmul(u, offset))
        at.pos = gemmi.Position(*p)
        res.add_atom(at)
    return res


def build(work: Path, prov: dict, numbering: pd.DataFrame, expected: str,
          jitter: dict[int, list[float]], drop_last: bool,
          break_at: int | None, drift_at: int | None) -> gemmi.Structure:
    chain_name = prov["chain"]
    src = gemmi.read_structure(str(work / "input.pdb"))
    src.setup_entities()
    obs = {i + 1: res for i, res in enumerate(src[0][chain_name])}

    rows = list(numbering.itertuples(index=False))
    # Pass 1: C-alpha for every full-pose position.
    ca: dict[int, list[float]] = {}
    for r in rows:
        if r.modelled == "0":
            atom = obs[int(r.input_pose)].find_atom("CA", "*")
            ca[int(r.full_pose)] = [atom.pos.x, atom.pos.y, atom.pos.z]

    truth_path = Path(prov["source_cif"]).with_suffix(".truth.tsv")
    if truth_path.exists():
        truth = pd.read_csv(truth_path, sep="\t", dtype=str)
        xyz = {row.auth_seq: [float(row.x), float(row.y), float(row.z)]
               for row in truth.itertuples(index=False)}
        for r in rows:
            if r.modelled == "1" and str(r.auth_seq) in xyz:
                ca[int(r.full_pose)] = xyz[str(r.auth_seq)]
        missing = [r for r in rows
                   if int(r.full_pose) not in ca]
        if missing:
            raise SystemExit(
                f"{truth_path} covers only part of the chain: "
                f"{len(missing)} position(s) have no C-alpha.")
        return finish(src, prov, rows, ca, jitter, drop_last, break_at,
                      drift_at)

    for g in prov["gaps"]:
        span = [int(r.full_pose) for r in rows
                if r.modelled == "1"
                and int(g["gap_start_auth"]) <= int(r.auth_seq) <= int(g["gap_end_auth"])]
        p0, p1 = ca[min(span) - 1], ca[max(span) + 1]
        normal = perpendicular(vnorm(vsub(p1, p0)))
        h = solve_bulge(p0, p1, len(span), normal)
        for pose, pt in zip(span, arc(p0, p1, len(span), h, normal)):
            ca[pose] = pt

    return finish(src, prov, rows, ca, jitter, drop_last, break_at, drift_at)


def finish(src, prov, rows, ca, jitter, drop_last, break_at, drift_at):
    chain_name = prov["chain"]
    for pose, off in jitter.items():
        if pose in ca:
            ca[pose] = vadd(ca[pose], off)
    if drift_at is not None:
        ca[drift_at] = vadd(ca[drift_at], [0.0, 5.0, 0.0])
    if break_at is not None:
        ca[break_at] = vadd(ca[break_at], [0.0, 12.0, 0.0])

    out = gemmi.Structure()
    out.name = "FAKE"
    out.cell = src.cell
    out.spacegroup_hm = src.spacegroup_hm or "P 1"
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_name)
    keep = rows[:-1] if drop_last else rows
    for r in keep:
        pose = int(r.full_pose)
        nxt = ca.get(pose + 1)
        direction = vsub(nxt, ca[pose]) if nxt else vsub(ca[pose], ca[pose - 1])
        chain.add_residue(residue(THREE[r.aa], pose, ca[pose], direction))
    model.add_chain(chain)
    for c in src[0]:
        if c.name != chain_name:
            copy = gemmi.Chain(c.name)
            for res in c:
                copy.add_residue(res)
            model.add_chain(copy)
    out.add_model(model)
    out.setup_entities()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--mode", choices=("insert", "ngk"), required=True)
    ap.add_argument("-n", type=int, default=3)
    ap.add_argument("--defect", choices=("none", "break", "sequence", "drift"),
                    default="none", help="plant this defect in candidate 1")
    ap.add_argument("--amplitude", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    work = Path(args.work_dir)
    prov = json.loads((work / "provenance.json").read_text())
    numbering = pd.read_csv(work / "numbering.tsv", sep="\t", dtype=str).fillna("")
    expected = (work / "expected_sequence.txt").read_text().strip()
    loop_poses = sorted({p for lp in prov["loops"]
                         for p in range(lp["lo"], lp["hi"] + 1)})
    interior = [p for p in loop_poses
                if p not in (lp["lo"] for lp in prov["loops"])
                and p not in (lp["hi"] for lp in prov["loops"])]
    fixed = [int(r.full_pose) for r in numbering.itertuples(index=False)
             if int(r.full_pose) not in loop_poses]
    rng = random.Random(args.seed)

    if args.mode == "insert":
        for k in range(1, args.n + 1):
            d = args.defect if k == 1 else "none"
            st = build(work, prov, numbering, expected, {},
                       drop_last=(d == "sequence"),
                       break_at=(interior[len(interior) // 2]
                                 if d == "break" else None),
                       drift_at=(fixed[0] if d == "drift" else None))
            st.write_pdb(str(work / f"{k}.pdb"))
        print(f"wrote {args.n} insert candidate(s), defect={args.defect}")
        return

    out = work / "models"
    out.mkdir(exist_ok=True)
    lines = ["SEQUENCE: ", "SCORE: total_score description"]
    for k in range(1, args.n + 1):
        jitter = {p: [rng.gauss(0, args.amplitude) for _ in range(3)]
                  for p in interior}
        name = f"ngk_{k:04d}"
        build(work, prov, numbering, expected, jitter,
              drop_last=False, break_at=None, drift_at=None
              ).write_pdb(str(out / f"{name}.pdb"))
        lines.append(f"SCORE: {-100.0 - k:.3f} {name}")
    # Beside the models, which is where check_models.py looks: each array
    # task writes its own score file next to its own output.
    (out / "ngk.sc").write_text("\n".join(lines) + "\n")
    print(f"wrote {args.n} ngk model(s), amplitude={args.amplitude}")


if __name__ == "__main__":
    main()
