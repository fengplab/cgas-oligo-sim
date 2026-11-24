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


def parse_header(header: str, seq_len: int | None = None,
                 coords: str = "0-half-open") -> dict:
    """Return the fields of one FASTA header. Raises HeaderError on anything odd.

    `coords` names the convention of the header's chrom:start-end token:

      "0-half-open"  end - start     == len(seq)   (BED/narrowPeak)
      "1-closed"     end - start + 1 == len(seq)   (samtools faidx region)

    Either way the returned start/end are **0-based half-open**, which is the
    only convention used downstream. A 1-closed header is converted here, at
    the single point where headers are read, so no other module has to know.

    crossrun's own R/plot_peak_sequences.R says the summit is 0-based
    (narrowPeak) while the faidx region it builds is 1-based closed, so a
    1-closed header carries a 0-based summit beside a 1-based interval. The
    summit is therefore NOT shifted, and the bounds check is done in 0-based
    space after conversion.
    """
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
    if coords == "1-closed":
        # faidx region -> 0-based half-open. Only the start moves.
        start, end = start - 1, end
    elif coords != "0-half-open":
        raise HeaderError(f"unknown coords convention {coords!r}")
    if not UCSC.match(chrom):
        raise HeaderError(f"'{oligo_id}': contig '{chrom}' is not UCSC-style. A "
                          "RefSeq (NC_0609xx.1) or GenBank (CP0682xx.2) name here "
                          "would fail to join the BAMs later.")
    if end <= start:
        raise HeaderError(f"'{oligo_id}': interval [{start},{end}) is empty or "
                          "reversed")
    if seq_len is not None and (end - start) != seq_len:
        alt = ("1-closed" if coords == "0-half-open" else "0-half-open")
        alt_span = (end - start + 1) if coords == "0-half-open" else (end - start - 1)
        hint = (f" Read as {alt} this record spans {alt_span} bp, which DOES "
                f"match; if that holds for every record in the file, pass "
                f"--header-coords {alt}."
                if alt_span == seq_len else
                f" Read as {alt} it would span {alt_span} bp, which does not "
                "match either, so this is not a convention difference.")
        raise HeaderError(
            f"'{oligo_id}': header interval spans {end - start} bp under "
            f"{coords} but the sequence is {seq_len} bp.{hint} Do not "
            "reconcile this by shifting one until you know which is wrong.")
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
    ap.add_argument("--diagnose-coords", action="store_true",
                    help="report, for every record, which coordinate "
                         "convention matches its sequence length, then exit "
                         "without writing anything.")
    ap.add_argument("--header-coords", choices=["0-half-open", "1-closed"],
                    default="0-half-open",
                    help="convention of the chrom:start-end token in the FASTA "
                         "headers. 0-half-open is BED/narrowPeak (end-start == "
                         "length); 1-closed is a samtools faidx region "
                         "(end-start+1 == length). crossrun's top-peaks FASTA "
                         "headers are faidx region strings, so they are "
                         "1-closed. Coordinates are stored 0-based half-open "
                         "either way.")
    ap.add_argument("--mutants", help="optional TSV: name, oligo, mutation")
    ap.add_argument("--min-length", type=int, default=60)
    ap.add_argument("--max-length", type=int, default=1000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-metadata", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if args.diagnose_coords:
        # One record agreeing with a convention is a coincidence; the whole
        # file agreeing is evidence. Writes nothing.
        recs = list(SeqIO.parse(args.fasta, "fasta"))
        if not recs:
            die(f"{args.fasta} contains no sequences")
        tally = {"0-half-open": 0, "1-closed": 0, "neither": 0}
        print(f"{'record':16s} {'len':>5s} {'end-start':>10s} {'+1':>5s}  fits")
        for rec in recs:
            seq = str(rec.seq)
            parts = (rec.description or rec.id).split()
            m = LOCUS.match(parts[1]) if len(parts) > 1 else None
            if not m:
                tally["neither"] += 1
                print(f"{parts[0]:16s} {len(seq):5d} {'-':>10s} {'-':>5s}  "
                      "no chrom:start-end token")
                continue
            half = int(m["end"]) - int(m["start"])
            fit = ("0-half-open" if half == len(seq) else
                   "1-closed" if half + 1 == len(seq) else "neither")
            tally[fit] += 1
            print(f"{parts[0]:16s} {len(seq):5d} {half:10d} {half+1:5d}  {fit}")
        n = len(recs)
        print()
        for k, v in tally.items():
            if v:
                print(f"  {v}/{n} record(s) fit {k}")
        fits = [k for k in ("0-half-open", "1-closed") if tally[k]]
        print()
        if tally["neither"]:
            print("Some records fit NEITHER convention, so this is not a "
                  "convention difference but a real inconsistency between the "
                  "headers and the sequences. Find it before building a panel "
                  "from this file.")
            sys.exit(1)
        if len(fits) == 1:
            print(f"Every record fits {fits[0]}. Rerun with "
                  f"--header-coords {fits[0]}.")
            sys.exit(0)
        print("The file MIXES both conventions. Do not pick one - the file is "
              "wrong and should be regenerated at source.")
        sys.exit(1)

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
            h = parse_header(rec.description or rec.id, seq_len=len(seq),
                             coords=args.header_coords)
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
            # What the input headers were READ as. Without this the stored
            # coordinates cannot be audited: a one-base shift leaves no trace
            # in the numbers themselves.
            "header_coords_input": args.header_coords,
            "locus": "0-based half-open, + strand (stored; converted from "
                     "header_coords_input)",
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
