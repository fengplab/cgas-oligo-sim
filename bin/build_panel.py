#!/usr/bin/env python3
"""Build the oligo panel from a coordinate-bearing FASTA and an optional mutant table.

FASTA header
------------
The oligo name is the first whitespace-delimited token; the rest of the
description carries the annotation:

  >rank1 chr1:130352560-130352639 summit=130352599 track=HSV-1 Infected THP1 Cells (Cytosolic Fraction - 1 HPI)

  name    rank1
  locus   chr1:130352560-130352639   0-based half-open, + strand
  summit  130352599                  optional, genomic
  track   everything after 'track=', free text, may contain spaces and parens

All sequences are taken to be on the + strand, so no reverse-complement
handling is applied anywhere.

Coordinate conventions, and the one place they differ
-----------------------------------------------------
Locus coordinates are **0-based half-open**, so end - start must equal the
sequence length; that is asserted per record.

Mutation positions are **1-based** into the oligo, because "the 30th base
pair" is an ordinal. These two conventions coexist in one panel, which is
where off-by-one errors breed, so every mutation's reference base is checked
against the parent sequence. A mismatch fails and reports whether the
0-based interpretation would have matched instead.

Mutant table (optional)
-----------------------
TSV with three columns:

  name       label for the mutant in results
  oligo      name of the parent sequence it derives from
  mutation   e.g. "30T>G", or several separated by ';': "30T>G;50A>C"

Without --mutants, no mutant is generated and the panel is the FASTA alone.
Mutants carry no coordinates of their own: they inherit the parent's through
the lineage, which is the only correct placement for a sequence that no
longer matches the genome.
"""
from __future__ import annotations

import argparse
import json
import re
import sys

import pandas as pd
from Bio import SeqIO

LOCUS = re.compile(r"(?P<chrom>chr[0-9XYM]+|[A-Za-z0-9._]+):(?P<start>\d+)-(?P<end>\d+)")
SUMMIT = re.compile(r"\bsummit=(\d+)")
TRACK = re.compile(r"\btrack=(.*)$")
MUT = re.compile(r"^(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
BASES = set("ACGT")


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def parse_header(desc: str, name: str) -> dict:
    """Pull locus, summit and track out of a FASTA description line."""
    rest = desc[len(name):].strip() if desc.startswith(name) else desc
    m = LOCUS.search(rest)
    if not m:
        die(f"{name}: no locus found in the header. Expected something like "
            f"'chr1:130352560-130352639' after the name. Got: {rest!r}")
    s, e = int(m.group("start")), int(m.group("end"))
    if e <= s:
        die(f"{name}: locus {m.group('chrom')}:{s}-{e} is empty or inverted")
    summit = SUMMIT.search(rest)
    track = TRACK.search(rest)
    return {"chrom": m.group("chrom"), "start": s, "end": e,
            "summit": int(summit.group(1)) if summit else None,
            "track": track.group(1).strip() if track else ""}


def apply_mutations(parent_seq: str, spec: str, mut_name: str,
                    parent_name: str, start: int) -> tuple[str, list[dict]]:
    """Apply 'POSREF>ALT[;...]' to parent_seq. Positions are 1-based."""
    seq = list(parent_seq)
    applied, touched = [], {}
    for piece in [p.strip() for p in spec.split(";") if p.strip()]:
        m = MUT.match(piece)
        if not m:
            die(f"{mut_name}: cannot parse mutation {piece!r}. Expected "
                "'<pos><REF>><ALT>' such as '30T>G', several separated by ';'.")
        pos, ref, alt = int(m.group("pos")), m.group("ref"), m.group("alt")
        if not 1 <= pos <= len(seq):
            die(f"{mut_name}: position {pos} is outside parent {parent_name}, "
                f"which is {len(seq)} bp (positions are 1-based)")
        if pos in touched:
            die(f"{mut_name}: position {pos} is mutated twice "
                f"({touched[pos]} and {piece}); the result would depend on order")
        actual = seq[pos - 1]
        if actual != ref:
            hint = ""
            if pos < len(seq) and seq[pos] == ref:
                hint = (f" Position {pos} is '{actual}', but {pos + 1} is '{ref}': "
                        "this looks like a 0-based position in a 1-based field.")
            die(f"{mut_name}: {piece} says the reference base at position {pos} "
                f"of {parent_name} is '{ref}', but it is '{actual}'.{hint}")
        if ref == alt:
            die(f"{mut_name}: {piece} changes '{ref}' to itself")
        seq[pos - 1] = alt
        touched[pos] = piece
        applied.append({"pos_1based": pos, "ref": ref, "alt": alt,
                        # + strand, 0-based half-open locus start
                        "genomic_0based": start + pos - 1})
    return "".join(seq), applied


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--mutants", help="optional TSV: name, oligo, mutation")
    ap.add_argument("--min-length", type=int, default=60)
    ap.add_argument("--max-length", type=int, default=1000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-metadata", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    parents, records, meta = {}, [], []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        name = rec.id
        seq = str(rec.seq).upper()
        if name in parents:
            die(f"duplicate oligo name in the FASTA: {name}")
        bad = sorted(set(seq) - BASES)
        if bad:
            die(f"{name} contains non-ACGT character(s) {bad}. Degenerate or "
                "masked bases cannot be built into a CG structure.")
        if not args.min_length <= len(seq) <= args.max_length:
            die(f"{name} is {len(seq)} bp, outside "
                f"[{args.min_length}, {args.max_length}]")
        h = parse_header(rec.description, name)
        width = h["end"] - h["start"]
        if width != len(seq):
            die(f"{name}: locus {h['chrom']}:{h['start']}-{h['end']} spans {width} bp "
                f"but the sequence is {len(seq)} bp. Coordinates are 0-based "
                "half-open; do not reconcile this by shifting one until you know "
                "which is wrong.")
        if h["summit"] is not None and not h["start"] <= h["summit"] < h["end"]:
            die(f"{name}: summit {h['summit']} lies outside the locus "
                f"[{h['start']}, {h['end']})")
        parents[name] = {"seq": seq, **h}
        records.append((name, seq))
        meta.append({
            "sequence_id": name, "lineage_id": name, "role": "wt",
            "chrom": h["chrom"], "start": h["start"], "end": h["end"], "strand": "+",
            "mutation": "",
            "notes": (f"summit={h['summit']}; " if h["summit"] is not None else "")
                     + (f"track={h['track']}" if h["track"] else ""),
        })

    if not parents:
        die(f"{args.fasta} contains no sequences")

    n_mut, mut_detail = 0, []
    if args.mutants:
        mt = pd.read_csv(args.mutants, sep="\t", dtype=str)
        need = ["name", "oligo", "mutation"]
        missing = [c for c in need if c not in mt.columns]
        if missing:
            die(f"{args.mutants} is missing column(s): {', '.join(missing)}. "
                f"Expected exactly {need}.")
        if mt["name"].duplicated().any():
            dup = sorted(mt.loc[mt["name"].duplicated(), "name"])[:5]
            die(f"duplicate mutant name(s): {dup}")
        clash = sorted(set(mt["name"]) & set(parents))
        if clash:
            die(f"mutant name(s) collide with a parent oligo: {clash[:5]}")

        for _, r in mt.iterrows():
            parent = r["oligo"]
            if parent not in parents:
                die(f"mutant {r['name']}: parent oligo '{parent}' is not in "
                    f"{args.fasta}. Available: "
                    + ", ".join(sorted(parents)[:8])
                    + (" ..." if len(parents) > 8 else ""))
            p = parents[parent]
            seq, applied = apply_mutations(p["seq"], str(r["mutation"]), r["name"],
                                           parent, p["start"])
            records.append((r["name"], seq))
            meta.append({
                "sequence_id": r["name"], "lineage_id": parent, "role": "mutant",
                "chrom": "", "start": "", "end": "", "strand": "",
                "mutation": r["mutation"],
                "notes": "inherits coordinates from parent " + parent,
            })
            n_mut += 1
            mut_detail.append({"name": r["name"], "parent": parent,
                               "changes": applied})

    with open(args.out_fasta, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pd.DataFrame(meta).to_csv(args.out_metadata, sep="\t", index=False)

    summary = {
        "n_parents": len(parents), "n_mutants": n_mut,
        "n_total": len(records),
        "mutants_supplied": bool(args.mutants),
        "length_min": min(len(s) for _, s in records),
        "length_max": max(len(s) for _, s in records),
        "conventions": {
            "locus": "0-based half-open, + strand",
            "mutation_position": "1-based into the oligo; reference base asserted",
        },
        "mutations": mut_detail,
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[build_panel] parents={len(parents)} mutants={n_mut} "
          f"total={len(records)} lengths={summary['length_min']}-"
          f"{summary['length_max']}bp")
    if not args.mutants:
        print("[build_panel] no --mutants supplied; mutation step skipped")


if __name__ == "__main__":
    main()
