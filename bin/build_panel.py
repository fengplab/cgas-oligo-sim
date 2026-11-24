#!/usr/bin/env python3
"""Build the oligo panel from a coordinate-bearing FASTA and an optional mutant table.

FASTA header
------------
This file holds the only header parser in the repository. Do not add a second
one anywhere: two parsers disagree eventually, and the disagreement surfaces
as a coordinate that is wrong rather than as an error.

The oligo name is the first whitespace-delimited token; the rest of the
description carries the annotation:

  >rank1 chr1:130352560-130352639 summit=130352599 track=HSV-1 Infected THP1 Cells (Cytosolic Fraction - 1 HPI)

A sequence that is deliberately not genomic - a published control oligo or a
designed construct - puts the literal token 'synthetic' where the locus goes:

  >G3_YSD synthetic role=control source=Herzner_2015_ni.3267_oligo122

Such records carry no coordinates, because inventing one would place a false
position in the metadata and in every later join. summit= is rejected for
them.

  name    rank1
  locus   chr1:130352560-130352639   0-based half-open, + strand
  summit  130352599                  optional, genomic
  track   everything after 'track=', free text, may contain spaces and parens

Also honoured, per parse_fasta_header.py: lineage=, role=, mutation=. Any
other key=value is kept verbatim in notes. Contigs must be UCSC-style; a
RefSeq or GenBank name is rejected, because it would fail to join the BAMs.

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

MUT = re.compile(r"^(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
BASES = set("ACGT")


# A key's value runs to the next " key=" or to end of line, so free-text
# values containing spaces (track=...) survive.
KV = re.compile(r"(\w+)=(.*?)(?=\s+\w+=|$)")
LOCUS = re.compile(r"^(?P<chrom>[^\s:]+):(?P<start>\d+)-(?P<end>\d+)$")
UCSC = re.compile(r"^chr([0-9]+|X|Y|M)$")


class HeaderError(ValueError):
    pass


def parse_header(header: str, seq_len: int | None = None) -> dict:
    """Return the fields of one FASTA header. Raises HeaderError on anything odd."""
    text = header[1:] if header.startswith(">") else header
    text = text.strip()
    if not text:
        raise HeaderError("empty header")

    parts = text.split(None, 2)
    oligo_id = parts[0]
    out = {"sequence_id": oligo_id, "lineage_id": oligo_id, "role": "wt",
           "chrom": None, "start": None, "end": None, "strand": "+",
           "summit": None, "mutation": "", "extra": {}}

    if len(parts) == 1:
        raise HeaderError(f"'{oligo_id}': header has no chrom:start-end token. "
                          "A sequence that is deliberately not genomic (a "
                          "published control oligo, a designed construct) uses "
                          "the literal token 'synthetic' in its place.")

    # A synthetic control oligo has no genomic locus, and inventing one would
    # put a false coordinate in the record and in any later join. Such records
    # carry no coordinates at all, exactly as mutants do.
    if parts[1] == "synthetic":
        out.update(chrom=None, start=None, end=None, role="synthetic")
        for key, val in KV.findall(parts[2] if len(parts) > 2 else ""):
            k, val = key.lower(), val.strip()
            if k == "summit":
                raise HeaderError(f"'{oligo_id}': summit= is meaningless for a "
                                  "synthetic sequence, which has no genomic "
                                  "position to be a summit of")
            if k == "lineage":
                out["lineage_id"] = val
            elif k == "role":
                out["role"] = val
            elif k == "mutation":
                out["mutation"] = val
            else:
                out["extra"][key] = val
        return out

    m = LOCUS.match(parts[1])
    if not m:
        raise HeaderError(f"'{oligo_id}': second token '{parts[1]}' is not "
                          "chrom:start-end")
    chrom, start, end = m["chrom"], int(m["start"]), int(m["end"])
    if not UCSC.match(chrom):
        raise HeaderError(f"'{oligo_id}': contig '{chrom}' is not UCSC-style. A "
                          "RefSeq (NC_0609xx.1) or GenBank (CP0682xx.2) name here "
                          "would fail to join the BAMs later.")
    if end <= start:
        raise HeaderError(f"'{oligo_id}': interval [{start},{end}) is empty or "
                          "reversed")
    if seq_len is not None and (end - start) != seq_len:
        raise HeaderError(
            f"'{oligo_id}': header interval [{start},{end}) spans {end - start} bp "
            f"but the sequence is {seq_len} bp. Coordinates are 0-based half-open; "
            "do not reconcile this by shifting one until you know which is wrong.")
    out.update(chrom=chrom, start=start, end=end)

    for key, val in KV.findall(parts[2] if len(parts) > 2 else ""):
        val = val.strip()
        k = key.lower()
        if k == "summit":
            try:
                s = int(val)
            except ValueError:
                raise HeaderError(f"'{oligo_id}': summit='{val}' is not an integer")
            if not start <= s < end:
                raise HeaderError(
                    f"'{oligo_id}': summit {s} lies outside the interval "
                    f"[{start},{end}). Either the summit or the interval is wrong, "
                    "and a summit-anchored window built from it would be silently "
                    "misplaced.")
            out["summit"] = s
        elif k == "lineage":
            out["lineage_id"] = val
        elif k == "role":
            out["role"] = val
        elif k == "mutation":
            out["mutation"] = val
        else:
            out["extra"][key] = val
    return out


def notes_from(parsed: dict) -> str:
    bits = []
    if parsed["summit"] is not None:
        bits.append(f"summit={parsed['summit']}")
    bits += [f"{k}={v}" for k, v in parsed["extra"].items()]
    return "; ".join(bits)


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


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
        seq = str(rec.seq).upper()
        bad = sorted(set(seq) - BASES)
        if bad:
            die(f"{rec.id} contains non-ACGT character(s) {bad}. Degenerate or "
                "masked bases cannot be built into a CG structure.")
        if not args.min_length <= len(seq) <= args.max_length:
            die(f"{rec.id} is {len(seq)} bp, outside "
                f"[{args.min_length}, {args.max_length}]")
        # Header parsing, the locus/length check, the UCSC contig assertion and
        # the summit bounds check all live in parse_fasta_header.
        try:
            h = parse_header(rec.description or rec.id, seq_len=len(seq))
        except HeaderError as exc:
            die(f"{args.fasta}: {exc}")
        name = h["sequence_id"]
        if name in parents:
            die(f"duplicate oligo name in the FASTA: {name}")
        parents[name] = {"seq": seq, **h}
        records.append((name, seq))
        meta.append({
            "sequence_id": name, "lineage_id": h["lineage_id"], "role": h["role"],
            "chrom": h["chrom"], "start": h["start"], "end": h["end"],
            "strand": h["strand"], "mutation": h["mutation"],
            "notes": notes_from(h),
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
