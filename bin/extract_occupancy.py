#!/usr/bin/env python3
"""Turn GENESIS CG trajectories into the occupancy tables the analysis consumes.

This is the step that defines what "occupied" means, so nothing downstream can
be more defensible than the choices here. All of them are parameters.

Definitions
-----------
contact        a cGAS bead within --contact-cutoff of any bead of nucleotide i
o(i)           per-nucleotide contact occupancy: the fraction of analysed
               frames in which nucleotide i is contacted by any cGAS copy.
               This is the UV-crosslink proxy, the motif site-occupancy
               quantity and the contact-map value - one definition, three
               uses, so they cannot disagree with each other.
binding event  a cGAS copy contacting a contiguous run of at least
               --footprint-min nucleotides. Shorter runs are glancing
               contacts, not binding, and are excluded from registers but
               still counted in o(i).
register       the midpoint of that contiguous run, 1-based
dwell_frac     per replicate, the fraction of analysed frames in which a
               register is occupied

Frames before --equil-frac of the trajectory are discarded.

Conventions
-----------
Nucleotide and register positions are **1-based indices into the oligo**,
matching the rest of the pipeline. They are not genomic coordinates.

Bead identity comes from --bead-map, never from position in the file. The map
must name, for every bead index in the trajectory, its molecule role
(dna|protein|partner), and for DNA its nucleotide index and strand, and for
protein its copy number.

A replicate in which cGAS never binds is a VALID result and is emitted with
zero occupancy everywhere. Only a malformed or unreadable trajectory fails.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def load_frames(path: str, n_beads_expected: int) -> np.ndarray:
    """Return coords (frames, beads, 3). DCD via MDAnalysis, or .npz for tests."""
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=False)
        if "coords" not in d:
            die(f"{path} has no 'coords' array")
        return np.asarray(d["coords"], dtype=np.float64)
    try:
        import MDAnalysis as mda
    except ImportError:
        die("MDAnalysis is required to read DCD trajectories")
    u = mda.Universe.empty(n_beads_expected, trajectory=True)
    u.load_new(path, format="DCD")
    return np.array([ts.positions.copy() for ts in u.trajectory], dtype=np.float64)


def contiguous_runs(idx: np.ndarray) -> list[tuple[int, int]]:
    if len(idx) == 0:
        return []
    runs, start, prev = [], idx[0], idx[0]
    for p in idx[1:]:
        if p == prev + 1:
            prev = p
            continue
        runs.append((start, prev))
        start, prev = p, p
    runs.append((start, prev))
    return runs


def principal_axis(pts: np.ndarray) -> np.ndarray:
    c = pts - pts.mean(axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return vt[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trajectory", required=True, help="DCD, or NPZ for testing")
    ap.add_argument("--bead-map", required=True,
                    help="TSV: bead_index, role(dna|protein|partner), nt_index, strand, copy")
    ap.add_argument("--sequence-id", required=True)
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--arm", default="baseline", help="arm label; never pool across arms")
    ap.add_argument("--contact-cutoff", type=float, default=10.0,
                    help="Angstrom. UNVERIFIED default: calibrate from the contact "
                         "distance distribution of a bound reference structure.")
    ap.add_argument("--footprint-min", type=int, default=12,
                    help="contiguous nucleotides contacted before it counts as binding")
    ap.add_argument("--equil-frac", type=float, default=0.2)
    ap.add_argument("--min-frames", type=int, default=50,
                    help="fail if fewer analysed frames than this")
    ap.add_argument("--out-occupancy", required=True)
    ap.add_argument("--out-contacts", required=True)
    ap.add_argument("--out-orientation", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if not 0 <= args.equil_frac < 1:
        die(f"--equil-frac must be in [0,1), got {args.equil_frac}")
    if args.contact_cutoff <= 0:
        die(f"--contact-cutoff must be > 0, got {args.contact_cutoff}")
    if args.footprint_min < 1:
        die(f"--footprint-min must be >= 1, got {args.footprint_min}")

    bm = pd.read_csv(args.bead_map, sep="\t")
    need = ["bead_index", "role"]
    missing = [c for c in need if c not in bm.columns]
    if missing:
        die(f"{args.bead_map} is missing column(s): {', '.join(missing)}")
    for c in ("nt_index", "strand", "copy"):
        if c not in bm.columns:
            bm[c] = pd.NA
    bm = bm.sort_values("bead_index").reset_index(drop=True)
    if list(bm["bead_index"]) != list(range(len(bm))):
        die(f"{args.bead_map} bead_index must be 0..N-1 with no gaps")

    dna = bm[bm["role"] == "dna"]
    pro = bm[bm["role"] == "protein"]
    if len(dna) == 0:
        die("bead map contains no dna beads")
    if len(pro) == 0:
        die("bead map contains no protein beads")
    if dna["nt_index"].isna().any():
        die("every dna bead needs an nt_index")
    nt_index = dna["nt_index"].astype(int).to_numpy()
    n_nt = int(nt_index.max())
    if int(nt_index.min()) != 1:
        die(f"nt_index must be 1-based and start at 1; minimum seen is {nt_index.min()}")

    copies = pro["copy"].fillna(1).astype(int)
    copy_ids = sorted(copies.unique())

    coords = load_frames(args.trajectory, len(bm))
    if coords.ndim != 3 or coords.shape[2] != 3:
        die(f"trajectory must give (frames, beads, 3); got {coords.shape}")
    if coords.shape[1] != len(bm):
        die(f"trajectory has {coords.shape[1]} beads but the bead map has "
            f"{len(bm)}. Identity must come from the map, so this is a hard "
            "mismatch, not something to reconcile by truncation.")
    if not np.isfinite(coords).all():
        die("trajectory contains non-finite coordinates")

    n_total = coords.shape[0]
    start = int(round(n_total * args.equil_frac))
    frames = coords[start:]
    if len(frames) < args.min_frames:
        die(f"{len(frames)} analysed frames after discarding {args.equil_frac:.0%} "
            f"equilibration (of {n_total} total) is below --min-frames "
            f"{args.min_frames}")

    from scipy.spatial import cKDTree

    dna_rows = dna["bead_index"].to_numpy()
    nt_of_bead = nt_index
    contact_count = np.zeros(n_nt + 1, dtype=np.int64)
    register_count: dict[int, int] = {}
    orient_rows = []
    n_events = 0
    n_frames_bound = 0      # frames with >=1 binding event, by any copy
    n_frames_contact = 0    # frames with >=1 contact, binding or glancing

    for f, xyz in enumerate(frames):
        dna_xyz = xyz[dna_rows]
        tree = cKDTree(dna_xyz)
        touched_any = set()
        frame_had_event = False
        for cid in copy_ids:
            rows = pro.loc[copies == cid, "bead_index"].to_numpy()
            hits = tree.query_ball_point(xyz[rows], r=args.contact_cutoff)
            nts = {int(nt_of_bead[j]) for sub in hits for j in sub}
            if not nts:
                continue
            touched_any |= nts
            runs = contiguous_runs(np.array(sorted(nts)))
            lo, hi = max(runs, key=lambda r: r[1] - r[0])
            if (hi - lo + 1) < args.footprint_min:
                continue
            reg = (lo + hi) // 2
            register_count[reg] = register_count.get(reg, 0) + 1
            n_events += 1
            frame_had_event = True
            # Orientation: cosine between the cGAS principal axis and the local
            # DNA axis at the register. Sign distinguishes the two 180-degree
            # flips, which is the quantity that corresponds to "strand".
            span = np.flatnonzero((nt_of_bead >= lo) & (nt_of_bead <= hi))
            if len(span) >= 4:
                dax = principal_axis(dna_xyz[span])
                pax = principal_axis(xyz[rows])
                orient_rows.append({
                    "sequence_id": args.sequence_id, "arm": args.arm,
                    "replicate": args.replicate, "frame": f + start,
                    "copy": cid, "register": reg,
                    "cos_axis_angle": float(np.dot(pax, dax)),
                })
        if touched_any:
            n_frames_contact += 1
        if frame_had_event:
            n_frames_bound += 1
        for i in touched_any:
            contact_count[i] += 1

    n_an = len(frames)
    occ = pd.DataFrame({
        "sequence_id": args.sequence_id, "arm": args.arm,
        "replicate": args.replicate,
        "position": np.arange(1, n_nt + 1),
        "dwell_frac": [register_count.get(p, 0) / n_an for p in range(1, n_nt + 1)],
    })
    con = pd.DataFrame({
        "sequence_id": args.sequence_id, "arm": args.arm,
        "replicate": args.replicate,
        "nt_position": np.arange(1, n_nt + 1),
        "contact_frac": contact_count[1:] / n_an,
    })
    occ.to_csv(args.out_occupancy, sep="\t", index=False)
    con.to_csv(args.out_contacts, sep="\t", index=False)
    pd.DataFrame(orient_rows, columns=["sequence_id", "arm", "replicate", "frame",
                                       "copy", "register", "cos_axis_angle"]
                 ).to_csv(args.out_orientation, sep="\t", index=False)

    summary = {
        "sequence_id": args.sequence_id, "arm": args.arm,
        "replicate": args.replicate,
        "frames_total": n_total, "frames_analysed": n_an,
        "equil_discarded": start,
        "n_nucleotides": n_nt, "n_protein_copies": len(copy_ids),
        "binding_events": n_events,
        # events_per_frame can exceed 1 when several cGAS copies bind at once,
        # so it is not a fraction and must not be read as one.
        "events_per_frame": float(n_events / n_an) if n_an else 0.0,
        "frames_with_binding_frac": float(n_frames_bound / n_an) if n_an else 0.0,
        "any_contact_frac": float(n_frames_contact / n_an) if n_an else 0.0,
        "no_binding_events": bool(n_events == 0),
        "contact_without_binding": bool(n_events == 0 and n_frames_contact > 0),
        "parameters": {"contact_cutoff": args.contact_cutoff,
                       "footprint_min": args.footprint_min,
                       "equil_frac": args.equil_frac},
        "conventions": "nucleotide and register positions are 1-based into the oligo",
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    if n_events:
        note = ""
    elif n_frames_contact:
        note = ("  NO BINDING EVENTS, contact only (valid): longest contiguous "
                f"contact never reached --footprint-min {args.footprint_min}")
    else:
        note = "  NEVER CONTACTED (valid)"
    print(f"[extract_occupancy] {args.sequence_id} rep{args.replicate} arm={args.arm} "
          f"frames={n_an}/{n_total} events={n_events} "
          f"bound_frames={summary['frames_with_binding_frac']:.3f} "
          f"contact_frames={summary['any_contact_frac']:.3f}" + note)


if __name__ == "__main__":
    main()
