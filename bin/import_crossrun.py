#!/usr/bin/env python3
"""Adapt fengplab/chipseq_crossrun outputs into this pipeline's inputs.

Mapping
-------
  crossrun work/seqs/<sample>.fa        -> oligo FASTA (one record per peak)
  crossrun work/seqs/<sample>.scan.bed  -> genomic coordinates (authoritative)
  crossrun manifest.tsv                 -> per-sample metadata carried onto each oligo
  crossrun work/motifs/family_reps.meme -> convert_motifs.py --meme  (used directly)

Why coordinates come from the BED and not the FASTA header
----------------------------------------------------------
The crossrun step-03 comment documents the header as
`<sample>|<peak>::chrom:start-end`, but `getfasta_canonical` in bin/lib.sh pipes
bedtools through `awk '/^>/ { sub(/::.*$/, "") }'`, which strips the coordinate
suffix. Depending on the bedtools version and that awk, the header may or may
not carry coordinates. `<sample>.scan.bed` always does, so it is the join key
here. Names are matched exactly between the two files and any mismatch is a
hard failure, never a silent drop.

Coordinates and contigs
-----------------------
scan.bed is BED: **0-based half-open**, which is the convention this pipeline's
metadata already uses, so no conversion is applied. crossrun normalises every
contig to UCSC style (`chr1`..`chrM`); that naming is asserted here so a RefSeq
or GenBank-named input cannot pass through unnoticed and collide with the BAMs
later.

Lineages
--------
Each imported peak becomes its own lineage with role `wt`, carrying the genomic
coordinates. Mutants are added afterwards as extra FASTA records sharing the
parent's `lineage_id` with role `mutant` and no coordinates of their own -
they inherit the parent's, which is the only way a mutant can be placed on the
genome, since by construction it does not match it.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

import pandas as pd
from Bio import SeqIO

UCSC = re.compile(r"^chr([0-9]+|X|Y|M)$")


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun-work", required=True,
                    help="crossrun <outdir>/work directory")
    ap.add_argument("--manifest", required=True, help="crossrun manifest.tsv")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="sample_id(s) to import; default every sample in the manifest")
    ap.add_argument("--min-length", type=int, default=60)
    ap.add_argument("--max-length", type=int, default=1000)
    ap.add_argument("--carry", nargs="*",
                    default=["series_id", "virus", "cell_type", "hpi",
                             "cgas_construct", "genome_build"],
                    help="manifest columns to carry into the notes field")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-metadata", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    man = pd.read_csv(args.manifest, sep="\t", dtype=str)
    if "sample_id" not in man.columns:
        die(f"{args.manifest} has no sample_id column")
    builds = sorted(set(man["genome_build"].dropna())) if "genome_build" in man else []
    if len(builds) > 1:
        die(f"manifest mixes genome builds {builds}. Every oligo must be on one "
            "build; split the import.")

    wanted = args.samples or list(man["sample_id"])
    unknown = [s for s in wanted if s not in set(man["sample_id"])]
    if unknown:
        die(f"sample(s) not in the manifest: {', '.join(unknown)}")

    seqdir = os.path.join(args.crossrun_work, "seqs")
    if not os.path.isdir(seqdir):
        die(f"{seqdir} does not exist. Point --crossrun-work at the crossrun "
            "<outdir>/work directory (it contains seqs/, motifs/, peaks/).")

    records, meta_rows = [], []
    counts = {"samples": 0, "peaks_in": 0, "too_short": 0, "too_long": 0,
              "non_acgt": 0, "kept": 0}
    skipped = []

    for sample in wanted:
        fa = os.path.join(seqdir, f"{sample}.fa")
        bed = os.path.join(seqdir, f"{sample}.scan.bed")
        for p in (fa, bed):
            if not os.path.exists(p):
                die(f"missing {p}. Run crossrun step 3 (steps/03_sequences.sh) for "
                    f"sample {sample} first.")
        counts["samples"] += 1

        coords = pd.read_csv(bed, sep="\t", header=None,
                             names=["chrom", "start", "end", "name"],
                             usecols=[0, 1, 2, 3], dtype={"chrom": str, "name": str})
        bad_contigs = sorted(set(coords.loc[~coords["chrom"].str.match(UCSC), "chrom"]))
        if bad_contigs:
            die(f"{bed} has non-UCSC contig name(s) {bad_contigs[:5]}. crossrun "
                "normalises to chr-style; a RefSeq (NC_0609xx.1) or GenBank "
                "(CP0682xx.2) name here would silently fail to join the BAMs.")
        by_name = coords.set_index("name")
        if by_name.index.duplicated().any():
            dup = sorted(by_name.index[by_name.index.duplicated()])[:5]
            die(f"{bed} has duplicate peak name(s) {dup}")

        seen = set()
        for rec in SeqIO.parse(fa, "fasta"):
            counts["peaks_in"] += 1
            name = rec.id
            seen.add(name)
            if name not in by_name.index:
                die(f"{fa}: record '{name}' has no row in {bed}. The FASTA and the "
                    "scan BED must name the same peaks; a mismatch means one of "
                    "them is stale.")
            c = by_name.loc[name]
            seq = str(rec.seq).upper()
            n = len(seq)
            if set(seq) - set("ACGT"):
                counts["non_acgt"] += 1
                skipped.append({"name": name, "reason": "non-ACGT base (soft-masked "
                                                        "or N in the reference)"})
                continue
            if n < args.min_length:
                counts["too_short"] += 1
                skipped.append({"name": name, "reason": f"{n} bp < min {args.min_length}"})
                continue
            if n > args.max_length:
                counts["too_long"] += 1
                skipped.append({"name": name, "reason": f"{n} bp > max {args.max_length}"})
                continue
            if (int(c["end"]) - int(c["start"])) != n:
                die(f"{name}: scan.bed interval spans "
                    f"{int(c['end']) - int(c['start'])} bp but the FASTA record is "
                    f"{n} bp. BED is 0-based half-open; do not reconcile this by "
                    "shifting a coordinate until you know which file is wrong.")

            sid = re.sub(r"[^A-Za-z0-9._-]", "_", name)
            mrow = man[man["sample_id"] == sample].iloc[0]
            notes = "; ".join(f"{k}={mrow[k]}" for k in args.carry
                              if k in man.columns and pd.notna(mrow[k]))
            records.append((sid, seq))
            meta_rows.append({
                "sequence_id": sid, "lineage_id": sid, "role": "wt",
                "chrom": c["chrom"], "start": int(c["start"]), "end": int(c["end"]),
                "strand": ".", "mutation": "",
                "notes": f"crossrun_sample={sample}; " + notes,
            })
            counts["kept"] += 1

        orphan = sorted(set(by_name.index) - seen)
        if orphan:
            die(f"{bed} names {len(orphan)} peak(s) absent from {fa}, e.g. "
                f"{orphan[:5]}. One of the two files is stale; rerun crossrun step 3.")

    if not records:
        die("no peaks survived import; nothing to write")

    ids = [r[0] for r in records]
    if len(set(ids)) != len(ids):
        dup = sorted({i for i in ids if ids.count(i) > 1})[:5]
        die(f"duplicate sequence_id across samples: {dup}. Peak names collide "
            "between samples; prefix them per sample before importing.")

    with open(args.out_fasta, "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pd.DataFrame(meta_rows).to_csv(args.out_metadata, sep="\t", index=False)

    summary = {"counts": counts, "genome_build": builds[0] if builds else None,
               "skipped": skipped[:200], "n_skipped": len(skipped),
               "coordinate_convention": "0-based half-open, taken from scan.bed",
               "contig_convention": "UCSC chr-style, asserted"}
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[import_crossrun] samples={counts['samples']} peaks_in={counts['peaks_in']} "
          f"kept={counts['kept']} skipped={len(skipped)} "
          f"(short={counts['too_short']} long={counts['too_long']} "
          f"non_acgt={counts['non_acgt']}) build={builds[0] if builds else 'unset'}")
    if skipped:
        print(f"[import_crossrun] {len(skipped)} peak(s) skipped; see "
              f"{args.out_summary} for the reason on each")


if __name__ == "__main__":
    main()
