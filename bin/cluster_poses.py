#!/usr/bin/env python3
"""Select a representative pose per site by gromos/Daura RMSD clustering.

Method
------
The convention in the MD structural literature: cluster pooled frames on
RMSD with the gromos (Daura) algorithm, take the most populated cluster, and
use a REAL frame as the representative rather than an averaged centroid (an
average of cartesian coordinates is not a physical structure).

The gromos algorithm yields a real frame directly:
  1. superpose all frames (Kabsch) on the interface beads;
  2. compute the pairwise RMSD matrix;
  3. the frame with the most neighbours within `cutoff` becomes a cluster
     centre; it and its neighbours form the cluster and leave the pool;
  4. repeat until the pool is empty.

Frames are pooled across replicates before clustering, which is why the
replicate of origin is carried through and reported: a cluster drawn from one
replicate is not the same evidence as one drawn from fifty.

Every representative is emitted with its cluster population. A pose seen in
3% of frames must never be presented like one seen in 80%.
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


def kabsch_superpose(coords: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Superpose each frame onto ref. coords (n, m, 3), ref (m, 3)."""
    out = np.empty_like(coords)
    ref_c = ref - ref.mean(axis=0)
    for i in range(coords.shape[0]):
        x = coords[i] - coords[i].mean(axis=0)
        u, _, vt = np.linalg.svd(x.T @ ref_c)
        d = np.sign(np.linalg.det(u @ vt))
        rot = u @ np.diag([1.0, 1.0, d]) @ vt
        out[i] = x @ rot
    return out


def pairwise_rmsd(coords: np.ndarray) -> np.ndarray:
    """Full pairwise RMSD matrix for superposed frames."""
    n, m, _ = coords.shape
    flat = coords.reshape(n, m * 3)
    sq = (flat ** 2).sum(axis=1)
    d2 = sq[:, None] + sq[None, :] - 2.0 * (flat @ flat.T)
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2 / m)


def gromos(rmsd: np.ndarray, cutoff: float) -> list[tuple[int, np.ndarray]]:
    """Daura/gromos clustering. Returns [(centre_index, member_indices), ...]."""
    n = rmsd.shape[0]
    remaining = np.ones(n, dtype=bool)
    clusters = []
    while remaining.any():
        idx = np.flatnonzero(remaining)
        sub = rmsd[np.ix_(idx, idx)] <= cutoff
        counts = sub.sum(axis=1)
        best = int(np.argmax(counts))
        members = idx[sub[best]]
        clusters.append((int(idx[best]), members))
        remaining[members] = False
    return clusters


def rmsd_diagnostics(rmsd: np.ndarray, n_frames: int) -> dict:
    """Off-diagonal RMSD distribution and a cutoff sweep, for auditing the cutoff."""
    off = rmsd[np.triu_indices_from(rmsd, k=1)]
    pct = {f"p{q}": float(np.percentile(off, q)) for q in (0, 5, 25, 50, 75, 95, 100)}
    sweep = []
    lo, hi = float(off.min()), float(off.max())
    for cut in np.linspace(lo, hi, 9):
        cl = gromos(rmsd, float(cut))
        top = max(len(m) for _, m in cl)
        sweep.append({"cutoff": round(float(cut), 3), "n_clusters": len(cl),
                      "top_population": top,
                      "top_population_frac": round(top / n_frames, 4)})
    return {"pairwise_rmsd_percentiles": pct, "cutoff_sweep": sweep}


def write_pdb(path: str, coords: np.ndarray, labels: np.ndarray,
              bfactors: np.ndarray, title: str) -> None:
    """Write CG beads as PDB pseudo-atoms. B-factor carries contact fraction."""
    with open(path, "w") as fh:
        fh.write(f"TITLE     {title[:60]}\n")
        fh.write("REMARK   1 Coarse-grained beads rendered as pseudo-atoms.\n")
        fh.write("REMARK   1 B-factor column = per-bead contact fraction (0-1 scaled x100).\n")
        for i, (lab, xyz, b) in enumerate(zip(labels, coords, bfactors), start=1):
            kind = str(lab).split("_")[0][:3].upper()
            chain = "D" if kind.startswith("DNA") else "P"
            resname = "DNA" if chain == "D" else "PRO"
            fh.write(
                f"ATOM  {i:5d}  CA  {resname:>3s} {chain}{i:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                f"{1.00:6.2f}{min(max(b, 0.0), 1.0) * 100:6.2f}          C\n"
            )
        fh.write("END\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", required=True)
    ap.add_argument("--site-id", required=True)
    ap.add_argument("--sequence-id", required=True)
    ap.add_argument("--contacts", required=True)
    ap.add_argument("--rmsd-cutoff", type=float, required=True)
    ap.add_argument("--min-frames", type=int, required=True)
    ap.add_argument("--min-top-population-frac", type=float, default=0.05,
                    help="fail if the most populated cluster holds less than this "
                         "fraction of frames: the cutoff is below the data's own "
                         "noise floor and the representative would be arbitrary")
    ap.add_argument("--report-only", action="store_true",
                    help="emit the RMSD distribution and a cutoff sweep, then exit "
                         "without selecting a representative. Run this FIRST to "
                         "choose --rmsd-cutoff from the data rather than guessing.")
    ap.add_argument("--out-pdb", required=True)
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if args.rmsd_cutoff <= 0:
        die(f"--rmsd-cutoff must be > 0, got {args.rmsd_cutoff}")

    data = np.load(args.npz, allow_pickle=False)
    for key in ("coords", "bead_labels", "replicate"):
        if key not in data:
            die(f"{args.npz} is missing array '{key}'")
    coords = np.asarray(data["coords"], dtype=np.float64)
    labels = np.asarray(data["bead_labels"])
    reps = np.asarray(data["replicate"])

    if coords.ndim != 3 or coords.shape[2] != 3:
        die(f"coords must have shape (frames, beads, 3), got {coords.shape}")
    n_frames = coords.shape[0]
    if n_frames < args.min_frames:
        die(f"site {args.site_id}: {n_frames} pooled frames is below --min-frames "
            f"({args.min_frames}); the representative pose would not be supported")
    if coords.shape[1] != len(labels):
        die(f"bead_labels length ({len(labels)}) does not match coords beads "
            f"({coords.shape[1]})")
    if len(reps) != n_frames:
        die(f"replicate array length ({len(reps)}) does not match frames ({n_frames})")
    if not np.isfinite(coords).all():
        die(f"site {args.site_id}: coords contain non-finite values")

    sup = kabsch_superpose(coords, coords[0])
    rmsd = pairwise_rmsd(sup)
    diag = rmsd_diagnostics(rmsd, n_frames)

    if args.report_only:
        with open(args.out_summary, "w") as fh:
            json.dump({"site_id": args.site_id, "sequence_id": args.sequence_id,
                       "mode": "report-only", "n_frames_pooled": n_frames,
                       **diag}, fh, indent=2)
        pct = diag["pairwise_rmsd_percentiles"]
        print(f"[cluster_poses] REPORT ONLY site={args.site_id} frames={n_frames}")
        print(f"  pairwise RMSD  min={pct['p0']:.2f}  p25={pct['p25']:.2f}  "
              f"median={pct['p50']:.2f}  p75={pct['p75']:.2f}  max={pct['p100']:.2f} A")
        print("  cutoff sweep (choose a cutoff where clusters are few AND the top "
              "population is a real majority):")
        for row in diag["cutoff_sweep"]:
            print(f"    {row['cutoff']:>7.2f} A -> {row['n_clusters']:>4d} clusters, "
                  f"top {row['top_population']:>4d} ({row['top_population_frac']:.0%})")
        print("  No representative written. Re-run with --rmsd-cutoff set.")
        return

    clusters = gromos(rmsd, args.rmsd_cutoff)

    # Degenerate regime 1: cutoff below the data's noise floor. Every frame is
    # its own cluster, so the "representative" is an arbitrary frame. This does
    # not error on its own, so it is caught here.
    top_frac = max(len(m) for _, m in clusters) / n_frames
    if top_frac < args.min_top_population_frac:
        pct = diag["pairwise_rmsd_percentiles"]
        die(f"site {args.site_id}: --rmsd-cutoff {args.rmsd_cutoff} A yields "
            f"{len(clusters)} clusters with the largest holding only {top_frac:.1%} "
            f"of {n_frames} frames (floor {args.min_top_population_frac:.0%}). The "
            f"cutoff is below this data's noise floor: minimum pairwise RMSD is "
            f"{pct['p0']:.2f} A, median {pct['p50']:.2f} A. Re-run with "
            f"--report-only to pick a cutoff from the data.")

    rows = []
    for rank, (centre, members) in enumerate(clusters, start=1):
        rows.append({
            "sequence_id": args.sequence_id,
            "site_id": args.site_id,
            "cluster_rank": rank,
            "centre_frame": centre,
            "population": len(members),
            "population_frac": len(members) / n_frames,
            "n_replicates_in_cluster": int(pd.unique(reps[members]).size),
            "mean_rmsd_to_centre": float(rmsd[centre, members].mean()),
        })
    cl = pd.DataFrame(rows)
    cl.to_csv(args.out_clusters, sep="\t", index=False)

    top = cl.iloc[0]
    centre_frame = int(top["centre_frame"])

    # Per-bead contact fraction -> B-factor. DNA beads are ordered along the
    # footprint, so bead k maps to the k-th position in the site's contact window.
    ct = pd.read_csv(args.contacts, sep="\t")
    sel = ct[(ct["site_id"] == args.site_id)]
    if len(sel) == 0:
        die(f"no contact rows for site {args.site_id} in {args.contacts}; the "
            "join produced zero overlap")
    per_pos = sel.groupby("position")["contact_frac"].max().sort_index()
    dna_idx = [i for i, l in enumerate(labels) if str(l).startswith("DNA")]
    if len(dna_idx) == 0:
        die("no bead labels start with 'DNA'; cannot map contact fractions")
    bf = np.zeros(len(labels))
    vals = per_pos.to_numpy()
    for j, bead in enumerate(dna_idx):
        if j < len(vals):
            bf[bead] = float(vals[j])

    write_pdb(args.out_pdb, sup[centre_frame], labels, bf,
              f"{args.site_id} gromos centre frame {centre_frame}")

    summary = {
        "site_id": args.site_id,
        "sequence_id": args.sequence_id,
        "method": "gromos (Daura) RMSD clustering; representative = cluster centre frame",
        "rmsd_cutoff_angstrom": args.rmsd_cutoff,
        "n_frames_pooled": n_frames,
        "n_replicates_pooled": int(pd.unique(reps).size),
        "n_clusters": len(clusters),
        "top_cluster_population": int(top["population"]),
        "top_cluster_population_frac": float(top["population_frac"]),
        "top_cluster_n_replicates": int(top["n_replicates_in_cluster"]),
        "representative_frame": centre_frame,
        "bfactor_meaning": "per-bead contact fraction x100 (DNA beads only; protein beads 0)",
        # Degenerate regime 2: one cluster holding every frame. This can be
        # genuine convergence, so it is flagged rather than failed - but a
        # 100% population claim must never be read as structural agreement
        # without checking the cutoff against the RMSD distribution below.
        "degenerate_single_cluster": bool(len(clusters) == 1
                                          and float(top["population_frac"]) >= 1.0),
        **diag,
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[cluster_poses] site={args.site_id} frames={n_frames} "
          f"clusters={len(clusters)} top_pop={int(top['population'])} "
          f"({top['population_frac']:.1%}) reps_in_top={int(top['n_replicates_in_cluster'])} "
          f"cutoff={args.rmsd_cutoff}")
    if summary["degenerate_single_cluster"]:
        print(f"[cluster_poses] WARNING site={args.site_id}: a single cluster holds "
              f"all {n_frames} frames at cutoff {args.rmsd_cutoff} A. This may be "
              f"genuine convergence, or a cutoff that has absorbed all variation. "
              f"Check pairwise_rmsd_percentiles in the summary before using it.")


if __name__ == "__main__":
    main()
