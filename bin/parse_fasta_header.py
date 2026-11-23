#!/usr/bin/env python3
"""Parse self-describing oligo FASTA headers.

Expected form
-------------
    >rank1 chr1:130352560-130352639 summit=130352599 track=HSV-1 Infected THP1 ...

  token 1    the oligo id. This becomes `sequence_id` (here: rank1).
  token 2    chrom:start-end in the reference.
  then       zero or more key=value pairs. A value may contain spaces: it runs
             until the next ` key=` or the end of the line, which is what lets
             `track=` carry a free-text description.

Conventions, fixed by the project
---------------------------------
  * Coordinates are **0-based half-open**, as BED is. `end - start` must equal
    the sequence length, and that is checked rather than assumed.
  * Sequences are always the **+ strand**, so strand is recorded as `+` and no
    reverse-complement handling exists anywhere.
  * Contigs are UCSC-style (`chr1`..`chrM`), matching crossrun's normalisation.

Recognised keys
---------------
  summit=   integer reference position; must fall inside [start, end)
  lineage=  lineage_id; defaults to the oligo id (each sequence its own lineage)
  role=     wt | mutant | control | shuffle; defaults to wt
  mutation= free text
  anything else is kept verbatim in `notes`.
"""
from __future__ import annotations

import re

# A key's value runs to the next " key=" or to end of line, so free-text values
# containing spaces (track=...) survive. A value containing its own "word=" would
# be split; no current producer does that.
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
        raise HeaderError(f"'{oligo_id}': header has no chrom:start-end token")

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
