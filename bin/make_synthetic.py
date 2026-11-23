#!/usr/bin/env python3
"""Generate synthetic test fixtures for the visualisation branch.

This exists to exercise the pipeline's logic without real trajectories. It
plants a known site in one sequence and leaves a second sequence as pure
background, so both the positive and the null path are covered.

Synthetic data only. Nothing here is a biological claim.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

BEADS = 60          # interface beads tracked per pose
FOOTPRINT = 20      # bp, cGAS register footprint (for contact map extent)


def occupancy_rows(rng, seq, length, n_rep, planted, strength, rep_hit_frac, shuffle_id=None):
    rows = []
    for rep in range(1, n_rep + 1):
        # Background: each replicate transiently visits a handful of registers.
        d = rng.beta(0.4, 30.0, size=length)
        n_visit = rng.integers(2, 6)
        for p in rng.choice(length, size=n_visit, replace=False):
            d[p] += rng.beta(2.0, 8.0) * 0.5
        for site in (planted or []):
            if rng.random() < rep_hit_frac:
                lo = max(0, site - 2)
                hi = min(length, site + 3)
                d[lo:hi] += strength * rng.uniform(0.7, 1.0)
        d = np.clip(d, 0.0, 1.0)
        for p in range(length):
            r = {"sequence_id": seq, "replicate": rep, "position": p + 1,
                 "dwell_frac": round(float(d[p]), 6)}
            if shuffle_id is not None:
                r["shuffle_id"] = shuffle_id
            rows.append(r)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--length", type=int, default=150)
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--shuffles", type=int, default=20)
    ap.add_argument("--frames", type=int, default=300)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # (seq, planted 0-based indices, strength, fraction of replicates hitting)
    specs = [
        ("SEQ_A", [59, 104], 0.45, 0.85),  # two planted sites, positions 60 and 105
        ("SEQ_B", [], 0.0, 0.0),           # background only
    ]

    obs, nul = [], []
    for seq, planted, strength, hitfrac in specs:
        obs += occupancy_rows(rng, seq, args.length, args.replicates,
                              planted, strength, hitfrac)
        for s in range(1, args.shuffles + 1):
            nul += occupancy_rows(rng, seq, args.length, max(8, args.replicates // 8),
                                  [], 0.0, 0.0, shuffle_id=f"shuf{s:02d}")

    obs_df = pd.DataFrame(obs)
    nul_df = pd.DataFrame(nul)[["sequence_id", "shuffle_id", "replicate", "position", "dwell_frac"]]
    obs_df.to_csv(f"{args.outdir}/occupancy.tsv", sep="\t", index=False)
    nul_df.to_csv(f"{args.outdir}/null_occupancy.tsv", sep="\t", index=False)

    # Poses: one NPZ per planted site. coords (frames, beads, 3) plus provenance.
    # Two conformational populations so clustering has something to separate.
    sheet = []
    for seq, planted, strength, _ in specs:
      for site in (planted or []):
        site_id = f"{seq}_s{site + 1}"
        base = rng.normal(0, 8, size=(BEADS, 3))
        alt = base + rng.normal(0, 3, size=(BEADS, 3))
        n_major = int(args.frames * 0.7)
        coords = np.concatenate([
            base[None, :, :] + rng.normal(0, 0.8, size=(n_major, BEADS, 3)),
            alt[None, :, :] + rng.normal(0, 0.8, size=(args.frames - n_major, BEADS, 3)),
        ])
        coords = coords[rng.permutation(args.frames)]
        reps = rng.integers(1, args.replicates + 1, size=args.frames)
        labels = np.array(
            [f"DNA_{i+1}" for i in range(FOOTPRINT)]
            + [f"PRO_{i+1}" for i in range(BEADS - FOOTPRINT)]
        )
        path = os.path.abspath(f"{args.outdir}/{site_id}_poses.npz")
        np.savez_compressed(path, coords=coords.astype(np.float32),
                            bead_labels=labels, replicate=reps)
        sheet.append({"sequence_id": seq, "site_id": site_id, "npz_path": path})

    pd.DataFrame(sheet).to_csv(f"{args.outdir}/poses.tsv", sep="\t", index=False)

    # Contact table: residue x position contact fractions around each site.
    crows = []
    for row in sheet:
        planted = int(row["site_id"].rsplit("_s", 1)[1])
        for res in range(1, 41):
            for off in range(-FOOTPRINT // 2, FOOTPRINT // 2 + 1):
                pos = planted + off
                if not 1 <= pos <= args.length:
                    continue
                # Zn-thumb proxy: residues 28-34 contact a fixed register offset.
                v = rng.beta(0.5, 12)
                if 28 <= res <= 34 and -3 <= off <= 1:
                    v += rng.beta(5, 3) * 0.6
                elif 5 <= res <= 20 and abs(off) <= 7:
                    v += rng.beta(2, 6) * 0.35
                crows.append({"sequence_id": row["sequence_id"], "site_id": row["site_id"],
                              "residue": res, "position": pos,
                              "contact_frac": round(float(min(v, 1.0)), 6)})
    pd.DataFrame(crows).to_csv(f"{args.outdir}/contacts.tsv", sep="\t", index=False)

    print(f"[make_synthetic] wrote fixtures to {args.outdir}")
    print(f"  occupancy.tsv       {len(obs_df):>8} rows")
    print(f"  null_occupancy.tsv  {len(nul_df):>8} rows  ({args.shuffles} shuffles/sequence)")
    print(f"  poses.tsv           {len(sheet):>8} sites")
    print(f"  contacts.tsv        {len(crows):>8} rows")


if __name__ == "__main__":
    main()
