#!/usr/bin/env python3
"""Call cGAS occupancy sites against a matched shuffle null.

Coordinate conventions
----------------------
`position` on input and in all TSV output is a **1-based index into the oligo
sequence**, naming the anchor of the cGAS register (the centre of the footprint,
as emitted by the upstream occupancy step). It is NOT a genomic coordinate.

`bed_start` / `bed_end` in sites.tsv are **0-based half-open** (BED convention),
derived as bed_start = position - 1, bed_end = position. They exist so the
sites table can be intersected with ChIP output without a conversion happening
somewhere implicit. Anything crossing into genomic space must first be lifted
onto the project's genome build; that is a separate step and is not done here.

Emptiness policy
----------------
Zero called sites is a scientifically valid outcome (it is the null result) and
is reported, not failed. Zero *input* rows, a sequence with no matched shuffles,
or an all-zero occupancy profile are failures, because each means the upstream
step did not produce what it claims to have produced.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

OCC_COLS = ["sequence_id", "replicate", "position", "dwell_frac"]
NULL_COLS = ["sequence_id", "shuffle_id", "replicate", "position", "dwell_frac"]
TIERS = (0.2, 0.4, 0.6, 0.8)


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def read_table(path: str, required: list[str], what: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in required if c not in df.columns]
    if missing:
        die(f"{what} ({path}) is missing required column(s): {', '.join(missing)}")
    if len(df) == 0:
        die(f"{what} ({path}) has no data rows")
    if df["dwell_frac"].isna().any():
        die(f"{what} ({path}) contains NaN in dwell_frac")
    bad = df.loc[(df["dwell_frac"] < 0) | (df["dwell_frac"] > 1), "dwell_frac"]
    if len(bad):
        die(f"{what} ({path}) has {len(bad)} dwell_frac value(s) outside [0,1]")
    if (df["position"] < 1).any():
        die(f"{what} ({path}) has position < 1; positions are 1-based")
    return df


def profile(df: pd.DataFrame, dwell_min: float) -> pd.DataFrame:
    """Per-position occupancy and replicate fraction for one sequence."""
    n_rep = df["replicate"].nunique()
    grouped = df.groupby("position", sort=True)
    out = pd.DataFrame(
        {
            "occupancy": grouped["dwell_frac"].mean(),
            "n_replicates_seen": grouped["replicate"].nunique(),
            "n_occupying": grouped["dwell_frac"].apply(lambda s: int((s >= dwell_min).sum())),
        }
    ).reset_index()
    out["replicate_frac"] = out["n_occupying"] / float(n_rep)
    return out


def bootstrap_ci(values: np.ndarray, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    """Percentile bootstrap CI of the mean, resampling replicates."""
    if len(values) < 2:
        return (float("nan"), float("nan"))
    draws = rng.choice(values, size=(n_boot, len(values)), replace=True).mean(axis=1)
    return (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)))


def merge_runs(positions: np.ndarray) -> list[tuple[int, int]]:
    """Collapse sorted called positions into contiguous [start, end] runs."""
    if len(positions) == 0:
        return []
    runs, start, prev = [], positions[0], positions[0]
    for p in positions[1:]:
        if p == prev + 1:
            prev = p
            continue
        runs.append((start, prev))
        start, prev = p, p
    runs.append((start, prev))
    return runs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--occupancy", required=True)
    ap.add_argument("--null-occupancy", required=True)
    ap.add_argument("--null-percentile", type=float, required=True)
    ap.add_argument("--dwell-min", type=float, required=True)
    ap.add_argument("--bootstrap-n", type=int, default=2000)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--replicate-min", default="", help="optional; empty means no replicate gate")
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--out-profile", required=True)
    ap.add_argument("--out-tiers", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if not 0 < args.null_percentile < 100:
        die(f"--null-percentile must be in (0,100), got {args.null_percentile}")
    if not 0 < args.dwell_min <= 1:
        die(f"--dwell-min must be in (0,1], got {args.dwell_min}")

    rep_min: float | None = None
    if str(args.replicate_min).strip() not in ("", "null", "None"):
        rep_min = float(args.replicate_min)
        if not 0 < rep_min <= 1:
            die(f"--replicate-min must be in (0,1] when set, got {rep_min}")

    rng = np.random.default_rng(args.seed)

    obs = read_table(args.occupancy, OCC_COLS, "occupancy table")
    nul = read_table(args.null_occupancy, NULL_COLS, "null occupancy table")

    obs_seqs = set(obs["sequence_id"].unique())
    nul_seqs = set(nul["sequence_id"].unique())
    unmatched = sorted(obs_seqs - nul_seqs)
    if unmatched:
        die("sequence(s) have no matched shuffle null, so no threshold can be set: "
            + ", ".join(map(str, unmatched[:10]))
            + (" ..." if len(unmatched) > 10 else ""))

    sites_rows, profile_rows, tier_rows = [], [], []
    counts = {"sequences_in": len(obs_seqs), "positions_in": 0, "positions_called": 0,
              "sites_out": 0, "sites_rendered": 0}

    for seq in sorted(obs_seqs):
        o = obs[obs["sequence_id"] == seq]
        n = nul[nul["sequence_id"] == seq]

        prof = profile(o, args.dwell_min)
        prof.insert(0, "sequence_id", seq)
        counts["positions_in"] += len(prof)

        if float(prof["occupancy"].max()) <= 0:
            die(f"sequence {seq}: occupancy is zero at every position. The upstream "
                "step produced no binding events; this is a failure, not 'no preference'.")

        # Null threshold: per-position occupancy computed within each shuffle,
        # then pooled across all positions and shuffles for this sequence.
        null_occ = (n.groupby(["shuffle_id", "position"], sort=False)["dwell_frac"]
                     .mean().to_numpy())
        if len(null_occ) == 0:
            die(f"sequence {seq}: matched null is empty after grouping")
        n_shuffles = int(n["shuffle_id"].nunique())
        threshold = float(np.percentile(null_occ, args.null_percentile))

        prof["null_threshold"] = threshold
        prof["above_threshold"] = prof["occupancy"] > threshold

        # Tier table: observed vs null counts at each replicate-fraction tier,
        # reported regardless of whether a replicate gate is in force. This is
        # what lets a replicate_min be chosen later from real data.
        null_prof = (n.groupby(["shuffle_id", "position"], sort=False)
                      .agg(occ=("dwell_frac", "mean"),
                           occupying=("dwell_frac", lambda s: int((s >= args.dwell_min).sum())),
                           nrep=("replicate", "nunique"))
                      .reset_index())
        null_prof["replicate_frac"] = null_prof["occupying"] / null_prof["nrep"].clip(lower=1)
        null_above = null_prof[null_prof["occ"] > threshold]

        for t in TIERS:
            tier_rows.append({
                "sequence_id": seq,
                "replicate_frac_tier": t,
                "observed_positions": int(((prof["occupancy"] > threshold)
                                           & (prof["replicate_frac"] >= t)).sum()),
                "null_positions_per_shuffle": (
                    float((null_above["replicate_frac"] >= t).sum()) / n_shuffles
                    if n_shuffles else float("nan")),
            })

        called = prof[prof["above_threshold"]]
        if rep_min is not None:
            called = called[called["replicate_frac"] >= rep_min]
        counts["positions_called"] += len(called)

        for start, end in merge_runs(called["position"].to_numpy()):
            window = prof[(prof["position"] >= start) & (prof["position"] <= end)]
            peak = window.loc[window["occupancy"].idxmax()]
            peak_pos = int(peak["position"])
            per_rep = o.loc[o["position"] == peak_pos, "dwell_frac"].to_numpy()
            lo, hi = bootstrap_ci(per_rep, args.bootstrap_n, rng)
            rendered = rep_min is None or float(peak["replicate_frac"]) >= rep_min
            sites_rows.append({
                "sequence_id": seq,
                "site_id": f"{seq}_s{peak_pos}",
                "peak_position": peak_pos,          # 1-based sequence index
                "bed_start": peak_pos - 1,          # 0-based half-open
                "bed_end": peak_pos,
                "run_start": int(start),
                "run_end": int(end),
                "width": int(end - start + 1),
                "occupancy": float(peak["occupancy"]),
                "occupancy_ci_lo": lo,
                "occupancy_ci_hi": hi,
                "replicate_frac": float(peak["replicate_frac"]),
                "n_occupying_replicates": int(peak["n_occupying"]),
                "null_threshold": threshold,
                "null_percentile": args.null_percentile,
                "n_shuffles": n_shuffles,
                "replicate_gate": "unset" if rep_min is None else f"{rep_min}",
                "rendered": bool(rendered),
            })

        profile_rows.append(prof)

    sites = pd.DataFrame(sites_rows)
    counts["sites_out"] = len(sites)
    counts["sites_rendered"] = int(sites["rendered"].sum()) if len(sites) else 0

    prof_all = pd.concat(profile_rows, ignore_index=True)
    prof_all.to_csv(args.out_profile, sep="\t", index=False)
    pd.DataFrame(tier_rows).to_csv(args.out_tiers, sep="\t", index=False)

    if len(sites) == 0:
        # Explicitly valid: write a header-only table and say so loudly.
        pd.DataFrame(columns=[
            "sequence_id", "site_id", "peak_position", "bed_start", "bed_end",
            "run_start", "run_end", "width", "occupancy", "occupancy_ci_lo",
            "occupancy_ci_hi", "replicate_frac", "n_occupying_replicates",
            "null_threshold", "null_percentile", "n_shuffles", "replicate_gate",
            "rendered"]).to_csv(args.out_sites, sep="\t", index=False)
    else:
        sites.sort_values(["sequence_id", "peak_position"]).to_csv(
            args.out_sites, sep="\t", index=False)

    summary = {
        "parameters": {
            "null_percentile": args.null_percentile,
            "dwell_min": args.dwell_min,
            "replicate_min": "unset" if rep_min is None else rep_min,
            "bootstrap_n": args.bootstrap_n,
            "seed": args.seed,
        },
        "counts": counts,
        "no_sites_called": bool(len(sites) == 0),
        "position_convention": "position is 1-based sequence index; bed_* are 0-based half-open",
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    # Counts in / counts out, to stdout for the run log.
    print(f"[call_sites] sequences={counts['sequences_in']} "
          f"positions_in={counts['positions_in']} "
          f"positions_called={counts['positions_called']} "
          f"sites_out={counts['sites_out']} "
          f"sites_rendered={counts['sites_rendered']} "
          f"replicate_gate={'unset' if rep_min is None else rep_min}")
    if len(sites) == 0:
        print("[call_sites] NO SITES CALLED. This is a valid null result: no position "
              f"exceeded the p{args.null_percentile:g} shuffle null.")


if __name__ == "__main__":
    main()
