#!/usr/bin/env python3
"""Synthetic CG trajectory with a KNOWN binding register, to test extract_occupancy.

Geometry: DNA laid along x as a straight ladder, 3 beads per nucleotide
(3SPN.2C-like), 3.4 A rise. One cGAS copy is a 60-bead blob that sits over a
chosen nucleotide window for a chosen fraction of frames and is parked far away
otherwise. The register the extractor reports must equal the window midpoint.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--n-nt", type=int, default=80)
ap.add_argument("--frames", type=int, default=200)
ap.add_argument("--bound-lo", type=int, default=28, help="1-based inclusive")
ap.add_argument("--bound-hi", type=int, default=45, help="1-based inclusive")
ap.add_argument("--bound-frac", type=float, default=0.6)
ap.add_argument("--copies", type=int, default=1)
ap.add_argument("--seed", type=int, default=20261123)
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
RISE, BPN = 3.4, 3
PRO_BEADS = 60

rows = []
bi = 0
for nt in range(1, a.n_nt + 1):
    for b in range(BPN):
        rows.append({"bead_index": bi, "role": "dna", "nt_index": nt,
                     "strand": 1 if b < 2 else 2, "copy": pd.NA})
        bi += 1
for c in range(1, a.copies + 1):
    for _ in range(PRO_BEADS):
        rows.append({"bead_index": bi, "role": "protein", "nt_index": pd.NA,
                     "strand": pd.NA, "copy": c})
        bi += 1
bm = pd.DataFrame(rows)
bm.to_csv(f"{a.out_prefix}.bead_map.tsv", sep="\t", index=False)

n_beads = len(bm)
dna_n = a.n_nt * BPN
dna_xyz = np.zeros((dna_n, 3))
for nt in range(a.n_nt):
    for b in range(BPN):
        k = nt * BPN + b
        dna_xyz[k] = [nt * RISE, 4.0 * (b - 1), 0.0]

coords = np.zeros((a.frames, n_beads, 3))
mid_x = ((a.bound_lo - 1) + (a.bound_hi - 1)) / 2.0 * RISE
half = ((a.bound_hi - a.bound_lo) / 2.0) * RISE
bound_flag = rng.random(a.frames) < a.bound_frac
for f in range(a.frames):
    coords[f, :dna_n] = dna_xyz + rng.normal(0, 0.3, dna_xyz.shape)
    for c in range(a.copies):
        s = dna_n + c * PRO_BEADS
        if bound_flag[f]:
            # slab hugging the window: within cutoff in y, spanning it in x
            xs = np.linspace(mid_x - half, mid_x + half, PRO_BEADS)
            blob = np.stack([xs, np.full(PRO_BEADS, 9.0), np.zeros(PRO_BEADS)], 1)
            blob += rng.normal(0, 0.4, blob.shape)
        else:
            blob = rng.normal(0, 3.0, (PRO_BEADS, 3)) + np.array([mid_x, 400.0, 0.0])
        coords[f, s:s + PRO_BEADS] = blob

np.savez_compressed(f"{a.out_prefix}.npz", coords=coords.astype(np.float32))
print(f"wrote {a.out_prefix}.npz  frames={a.frames} beads={n_beads} "
      f"expected_register={(a.bound_lo + a.bound_hi) // 2} "
      f"bound_frames={int(bound_flag.sum())}/{a.frames} "
      f"(post-equil truth depends on --equil-frac)")
