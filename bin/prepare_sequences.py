#!/usr/bin/env python3
"""Validate oligo FASTA + lineage metadata into a canonical samplesheet.

Why lineage metadata is load-bearing
------------------------------------
ChIP BAMs are aligned to the genome, but a mutant oligo differs from the genome
by construction, so it CANNOT be located by alignment. Genomic coordinates are
therefore declared once on each lineage's wild-type parent and inherited by its
mutants. A lineage with no WT member has no coordinates and is excluded from any
genome-anchored step, loudly.

Coordinates
-----------
Metadata `start`/`end` are **0-based half-open (BED)**, matching the BAM/BED
convention, and are checked against the sequence length. Positions emitted for
the oligo itself elsewhere in this pipeline are 1-based sequence indices; the
two are never mixed without an explicit conversion.

Oligo length is free to vary between sequences. Within one lineage, members are
expected to share a length (a substitution panel); an indel that changes length
is allowed but flagged, because it silently breaks positional comparison across
the lineage.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_fasta_header import HeaderError, notes_from, parse_header  # noqa: E402

META_REQUIRED = ["sequence_id", "lineage_id", "role"]
META_OPTIONAL = ["chrom", "start", "end", "strand", "mutation", "notes"]
VALID_ROLES = {"wt", "mutant", "control", "shuffle"}
IUPAC_DNA = set("ACGT")


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--metadata",
                    help="TSV: sequence_id, lineage_id, role [, chrom, start, end, "
                         "strand, mutation, notes]. Optional when --from-headers "
                         "is given, in which case it supplies only lineage_id, "
                         "role and mutation overrides.")
    ap.add_argument("--from-headers", action="store_true",
                    help="take sequence_id, coordinates and summit from the FASTA "
                         "headers: '>rank1 chr1:100-179 summit=140 track=...'. "
                         "Coordinates are 0-based half-open and + strand.")
    ap.add_argument("--min-length", type=int, default=60,
                    help="human cGAS responds preferentially to DNA >45 bp; a "
                         "shorter oligo sits where it is built not to respond")
    ap.add_argument("--max-length", type=int, default=1000)
    ap.add_argument("--require-coords", action="store_true",
                    help="fail if any lineage lacks genomic coordinates on its WT "
                         "member. Use when a genome-anchored step is downstream.")
    ap.add_argument("--out-samplesheet", required=True)
    ap.add_argument("--out-lineages", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if not args.from_headers and not args.metadata:
        die("give --metadata, or --from-headers to read coordinates out of the "
            "FASTA headers (or both: headers for coordinates, the TSV for "
            "lineage and role).")

    records, header_rows = {}, []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        seq = str(rec.seq).upper()
        if args.from_headers:
            raw = rec.description or rec.id
            try:
                h = parse_header(raw, seq_len=len(seq))
            except HeaderError as exc:
                die(f"{args.fasta}: {exc}")
            if h["sequence_id"] in records:
                die(f"duplicate oligo id '{h['sequence_id']}' in {args.fasta}")
            records[h["sequence_id"]] = seq
            header_rows.append({
                "sequence_id": h["sequence_id"], "lineage_id": h["lineage_id"],
                "role": h["role"], "chrom": h["chrom"], "start": str(h["start"]),
                "end": str(h["end"]), "strand": h["strand"],
                "mutation": h["mutation"], "notes": notes_from(h)})
        else:
            if rec.id in records:
                die(f"duplicate sequence id '{rec.id}' in {args.fasta}")
            records[rec.id] = seq
    if not records:
        die(f"{args.fasta} contains no sequences")

    if args.from_headers:
        meta = pd.DataFrame(header_rows).astype(str)
        if args.metadata:
            # The TSV supplies only what a header cannot: lineage grouping, role
            # and the mutation description. Coordinates always come from the
            # header, so the two can never disagree about position.
            over = pd.read_csv(args.metadata, sep="\t", dtype=str)
            if "sequence_id" not in over.columns:
                die(f"{args.metadata} has no sequence_id column")
            unknown = sorted(set(over["sequence_id"]) - set(meta["sequence_id"]))
            if unknown:
                die(f"{args.metadata} names sequence(s) absent from the FASTA: "
                    + ", ".join(unknown[:10]))
            meta = meta.set_index("sequence_id")
            for _, r in over.iterrows():
                for col in ("lineage_id", "role", "mutation"):
                    if col in over.columns and pd.notna(r.get(col)) and r.get(col) != "":
                        meta.loc[r["sequence_id"], col] = r[col]
            meta = meta.reset_index()
    else:
        meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    missing_cols = [c for c in META_REQUIRED if c not in meta.columns]
    if missing_cols:
        die(f"{args.metadata} is missing required column(s): {', '.join(missing_cols)}")
    for c in META_OPTIONAL:
        if c not in meta.columns:
            meta[c] = pd.NA

    if meta["sequence_id"].duplicated().any():
        dup = sorted(meta.loc[meta["sequence_id"].duplicated(), "sequence_id"])
        die(f"duplicate sequence_id in metadata: {', '.join(dup[:10])}")

    # Identity join: FASTA headers against the metadata, both directions.
    fa_ids, md_ids = set(records), set(meta["sequence_id"])
    only_fa, only_md = sorted(fa_ids - md_ids), sorted(md_ids - fa_ids)
    if only_fa:
        die(f"{len(only_fa)} sequence(s) in the FASTA have no metadata row: "
            + ", ".join(only_fa[:10]) + (" ..." if len(only_fa) > 10 else ""))
    if only_md:
        die(f"{len(only_md)} metadata row(s) name a sequence absent from the FASTA: "
            + ", ".join(only_md[:10]) + (" ..." if len(only_md) > 10 else ""))

    bad_roles = sorted(set(meta["role"]) - VALID_ROLES)
    if bad_roles:
        die(f"unknown role(s) {bad_roles}; valid roles are {sorted(VALID_ROLES)}")

    rows, warnings = [], []
    for _, m in meta.iterrows():
        sid = m["sequence_id"]
        seq = records[sid]
        bad = sorted(set(seq) - IUPAC_DNA)
        if bad:
            die(f"{sid} contains non-ACGT character(s) {bad}. Degenerate or masked "
                "bases cannot be built into a CG structure; resolve them upstream.")
        n = len(seq)
        if n < args.min_length:
            die(f"{sid} is {n} bp, below --min-length {args.min_length}. Human cGAS "
                "responds preferentially to DNA >45 bp, so a null here is uninformative.")
        if n > args.max_length:
            die(f"{sid} is {n} bp, above --max-length {args.max_length}")
        rows.append({
            "sequence_id": sid, "lineage_id": m["lineage_id"], "role": m["role"],
            "length": n, "sequence": seq,
            "chrom": m["chrom"], "start": m["start"], "end": m["end"],
            "strand": m["strand"] if pd.notna(m["strand"]) else ".",
            "mutation": m["mutation"], "notes": m["notes"],
        })

    df = pd.DataFrame(rows)

    # Coordinate sanity on whoever declares them.
    for _, r in df[df["chrom"].notna()].iterrows():
        if pd.isna(r["start"]) or pd.isna(r["end"]):
            die(f"{r['sequence_id']} declares chrom but not both start and end")
        s, e = int(r["start"]), int(r["end"])
        if s < 0 or e <= s:
            die(f"{r['sequence_id']} has invalid 0-based half-open interval "
                f"[{s},{e})")
        if (e - s) != r["length"]:
            die(f"{r['sequence_id']}: interval [{s},{e}) spans {e - s} bp but the "
                f"sequence is {r['length']} bp. BED is 0-based half-open; an "
                "off-by-one here silently shifts every downstream coordinate.")

    # Lineage roll-up. Coordinates are inherited from the WT parent.
    lin_rows = []
    for lid, g in df.groupby("lineage_id", sort=True):
        wt = g[g["role"] == "wt"]
        if len(wt) > 1:
            die(f"lineage {lid} has {len(wt)} members with role 'wt'; exactly one "
                "parent defines the lineage's coordinates")
        has_wt = len(wt) == 1

        # Coordinates for the lineage come from the wt member when it declares
        # them, otherwise from any member that does. Members that declare their
        # own keep them; the lineage value exists so that members declaring none
        # (a mutant written without a locus) can inherit.
        declared = g[g["chrom"].notna() & (g["chrom"].astype(str) != "")]
        loci = {(r["chrom"], int(r["start"]), int(r["end"]))
                for _, r in declared.iterrows()}
        if len(loci) > 1:
            die(f"lineage {lid} has members declaring different loci {sorted(loci)}. "
                "Members of one lineage occupy the same reference interval; a "
                "substitution mutant does not move. Fix the headers rather than "
                "letting one of them win.")
        coords = {}
        if loci:
            src = wt.iloc[0] if (has_wt and pd.notna(wt.iloc[0]["chrom"])
                                 and str(wt.iloc[0]["chrom"]) != "") \
                else declared.iloc[0]
            coords = {"chrom": src["chrom"], "start": int(src["start"]),
                      "end": int(src["end"]), "strand": src["strand"]}
        lengths = sorted(set(g["length"]))
        if len(lengths) > 1:
            warnings.append(
                f"lineage {lid} has members of differing length {lengths}: an indel "
                "breaks positional comparison across the lineage, so occupancy "
                "positions are not aligned between its members")
        if not has_wt:
            warnings.append(f"lineage {lid} has no 'wt' member, so it inherits no "
                            "genomic coordinates and is excluded from genome-anchored steps")
        elif not coords:
            warnings.append(f"lineage {lid} has a 'wt' member but it declares no "
                            "chrom/start/end, so the lineage is excluded from "
                            "genome-anchored steps (ChIP coverage, liftover)")
        lin_rows.append({
            "lineage_id": lid, "n_members": len(g), "has_wt": has_wt,
            "roles": ",".join(sorted(set(g["role"]))),
            "lengths": ",".join(map(str, lengths)),
            "genome_anchored": bool(coords),
            **{k: coords.get(k) for k in ("chrom", "start", "end", "strand")},
        })
    lineages = pd.DataFrame(lin_rows)
    for c in ("start", "end"):
        lineages[c] = lineages[c].astype("Int64")

    # Propagate inherited coordinates onto every member.
    inherit = lineages.set_index("lineage_id")[["chrom", "start", "end", "strand"]]
    inherit = inherit.astype({"start": "Int64", "end": "Int64"})
    df = df.join(inherit.rename(columns=lambda c: f"lineage_{c}"), on="lineage_id")
    # Effective position: what the record declares, else what it inherits.
    for c in ("chrom", "start", "end", "strand"):
        own = df[c].replace("", pd.NA)
        df[f"effective_{c}"] = own.where(own.notna(), df[f"lineage_{c}"])
    for c in ("start", "end"):
        df[f"effective_{c}"] = pd.to_numeric(df[f"effective_{c}"],
                                             errors="coerce").astype("Int64")

    if args.require_coords:
        un = sorted(lineages.loc[~lineages["genome_anchored"], "lineage_id"])
        if un:
            die("--require-coords is set but these lineage(s) have no WT genomic "
                "coordinates: " + ", ".join(un))

    df.to_csv(args.out_samplesheet, sep="\t", index=False)
    lineages.to_csv(args.out_lineages, sep="\t", index=False)
    summary = {
        "n_sequences": len(df),
        "n_lineages": len(lineages),
        "n_genome_anchored_lineages": int(lineages["genome_anchored"].sum()),
        "length_min": int(df["length"].min()),
        "length_max": int(df["length"].max()),
        "roles": df["role"].value_counts().to_dict(),
        "warnings": warnings,
        "coordinate_convention": "metadata start/end are 0-based half-open (BED)",
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[prepare_sequences] in={len(records)} out={len(df)} "
          f"lineages={len(lineages)} genome_anchored={summary['n_genome_anchored_lineages']} "
          f"lengths={summary['length_min']}-{summary['length_max']}bp")
    for w in warnings:
        print(f"[prepare_sequences] WARNING: {w}")


if __name__ == "__main__":
    main()
