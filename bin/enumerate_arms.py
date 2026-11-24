#!/usr/bin/env python3
"""Expand the arm configuration into one flat job matrix.

Every trajectory in the run is one row of jobs.tsv. Nextflow fans out over
the rows directly, so the whole run is a single wide queue with no nested
loops and no barrier between arms: the scheduler can keep as many jobs in
flight as the partition allows, for as long as the run takes.

The arms
--------
baseline      free association, no partner, no PWMcos restraint. Always present.
partner:<P>   free association with a DNA-binding partner (TFAM, HMGB1, ...).
              One arm per partner in --partners.
pwmcos:<M>    PWMcos sequence-recognition restraint for ONE motif family. A
              PWMcos term encodes one specific PWM, so arms cannot share it.
              **No motifs passed means no PWMcos arms at all** - the free
              association arms stand on their own and the run is far smaller.
structure:<S> each structure beyond the first gets its own single-basin arm,
              so a per-structure effect can be read without the multi-basin
              model in the way.
multibasin    one arm using all structures as basins, present only when two or
              more structures are given.

Observed rows use the real oligo panel. Null rows use shuffled sequences and
exist for every arm, because a threshold learned under one arm's force field
does not transfer to another's.

Job identity is built from the samplesheets, never from file order, and every
row carries an explicit seed so a single job can be reproduced on its own.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys

import pandas as pd


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def read_ids(path: str, cols: tuple[str, ...], what: str) -> list[str]:
    """IDs from the first of `cols` the file actually has.

    The panel is read from build_panel.py's metadata, whose key is
    `sequence_id` - the same key extract_occupancy.py and call_sites.py use, so
    one identifier runs the length of the pipeline. A bare one-column list
    headed `name` is also accepted for a hand-written panel.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    col = next((c for c in cols if c in df.columns), None)
    if col is None:
        die(f"{path} has none of the expected column(s) {list(cols)} "
            f"(needed to name {what}); it has: {list(df.columns)}")
    ids = df[col].dropna().astype(str).str.strip().tolist()
    if not ids:
        die(f"{path} lists no {what}")
    dups = sorted({i for i in ids if ids.count(i) > 1})
    if dups:
        die(f"{path} has duplicate {col} value(s): {dups}")
    return ids


def seed_for(global_seed: int, *parts: object) -> int:
    key = f"{global_seed}|" + "|".join(str(p) for p in parts)
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True,
                    help="build_panel.py's panel metadata TSV (key column "
                         "'sequence_id'); a bare list headed 'name' also works")
    ap.add_argument("--structures", required=True,
                    help="TSV of prepared structures; needs 'structure_id'")
    ap.add_argument("--motifs", default="",
                    help="optional TSV of motif families; needs 'family'. "
                         "Omit it and no PWMcos arm is generated.")
    ap.add_argument("--partners", default="",
                    help="comma-separated partner names, e.g. TFAM or "
                         "TFAM,HMGB1. Empty means no partner arm.")
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--shuffles", type=int, default=20)
    ap.add_argument("--null-replicates", type=int, default=12)
    ap.add_argument("--multibasin/--no-multibasin", dest="multibasin",
                    action="store_true", default=True)
    ap.add_argument("--no-multibasin", dest="multibasin", action="store_false")
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--out-jobs", required=True)
    ap.add_argument("--out-arms", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    for name, v in (("--replicates", args.replicates),
                    ("--shuffles", args.shuffles),
                    ("--null-replicates", args.null_replicates)):
        if v < 1:
            die(f"{name} must be >= 1, got {v}")

    oligos = read_ids(args.panel, ("sequence_id", "name"), "oligos")
    structures = read_ids(args.structures, ("structure_id",), "structures")
    motifs = read_ids(args.motifs, ("family",), "motif families") if args.motifs else []
    partners = [p.strip() for p in args.partners.split(",") if p.strip()]
    if len(set(partners)) != len(partners):
        die(f"--partners repeats a name: {partners}")

    arms = [{"arm": "baseline", "mode": "free", "partner": "",
             "motif": "", "structures": structures[0], "basins": 1}]
    for p in partners:
        arms.append({"arm": f"partner:{p}", "mode": "free", "partner": p,
                     "motif": "", "structures": structures[0], "basins": 1})
    for m in motifs:
        arms.append({"arm": f"pwmcos:{m}", "mode": "pwmcos", "partner": "",
                     "motif": m, "structures": structures[0], "basins": 1})
    for s in structures[1:]:
        arms.append({"arm": f"structure:{s}", "mode": "free", "partner": "",
                     "motif": "", "structures": s, "basins": 1})
    if len(structures) >= 2 and args.multibasin:
        arms.append({"arm": "multibasin", "mode": "free", "partner": "",
                     "motif": "", "structures": ",".join(structures),
                     "basins": len(structures)})

    arms_df = pd.DataFrame(arms)
    arms_df.to_csv(args.out_arms, sep="\t", index=False)

    rows = []
    for a in arms:
        for o in oligos:
            for r in range(1, args.replicates + 1):
                rows.append({
                    "job_id": f"{a['arm']}|{o}|obs|r{r}",
                    "arm": a["arm"], "mode": a["mode"], "partner": a["partner"],
                    "motif": a["motif"], "structures": a["structures"],
                    "basins": a["basins"], "sequence_id": o,
                    "panel": "observed", "shuffle": 0, "replicate": r,
                    "seed": seed_for(args.seed, a["arm"], o, "obs", r),
                })
            for sh in range(1, args.shuffles + 1):
                for r in range(1, args.null_replicates + 1):
                    rows.append({
                        "job_id": f"{a['arm']}|{o}|null{sh}|r{r}",
                        "arm": a["arm"], "mode": a["mode"],
                        "partner": a["partner"], "motif": a["motif"],
                        "structures": a["structures"], "basins": a["basins"],
                        "sequence_id": o, "panel": "null", "shuffle": sh,
                        "replicate": r,
                        "seed": seed_for(args.seed, a["arm"], o, f"null{sh}", r),
                    })

    jobs = pd.DataFrame(rows)
    if jobs["job_id"].duplicated().any():
        d = jobs.loc[jobs["job_id"].duplicated(), "job_id"].head().tolist()
        die(f"job_id collision, which would make jobs overwrite each "
            f"other's output: {d}")
    jobs.to_csv(args.out_jobs, sep="\t", index=False)

    n_obs = int((jobs.panel == "observed").sum())
    n_null = int((jobs.panel == "null").sum())
    summary = {
        "n_arms": len(arms),
        "arms": [a["arm"] for a in arms],
        "pwmcos_arms": len(motifs),
        "pwmcos_skipped_no_motifs": bool(not motifs),
        "multibasin_arm": bool(len(structures) >= 2 and args.multibasin),
        "n_oligos": len(oligos), "n_structures": len(structures),
        "n_partners": len(partners), "n_motif_families": len(motifs),
        "replicates": args.replicates, "shuffles": args.shuffles,
        "null_replicates": args.null_replicates,
        "jobs_observed": n_obs, "jobs_null": n_null,
        "jobs_total": len(jobs),
        "null_fraction": round(n_null / len(jobs), 4) if len(jobs) else 0.0,
        "seed": args.seed,
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[enumerate_arms] arms={len(arms)}  "
          f"({'no PWMcos arms: no motifs passed' if not motifs else f'{len(motifs)} PWMcos'}"
          f", {len(partners)} partner, {len(structures)} structure"
          f"{', multibasin' if summary['multibasin_arm'] else ''})")
    for a in arms:
        print(f"    {a['arm']:28s} mode={a['mode']:7s} basins={a['basins']}")
    print(f"[enumerate_arms] jobs: {n_obs:,} observed + {n_null:,} null "
          f"= {len(jobs):,} total  (null {summary['null_fraction']:.0%})")


if __name__ == "__main__":
    main()
