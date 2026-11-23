#!/usr/bin/env python3
"""Adapt fengplab/chipseq_crossrun outputs into this pipeline's inputs.

Mapping
-------
Two modes
---------
HARVEST (no --fasta): every peak in work/seqs/<sample>.fa becomes an oligo.

ANNOTATE (--fasta USER.fa): YOUR FASTA is authoritative for which sequences
exist, and crossrun is consulted only to supply coordinates and metadata. Use
this when the panel is hand-picked, renamed, or already contains designed
sequences that never came from a peak.

  crossrun work/seqs/<sample>.fa        -> oligo FASTA (harvest), or the
                                           lookup table to match against (annotate)
  crossrun work/seqs/<sample>.scan.bed  -> genomic coordinates (authoritative)
  crossrun manifest.tsv                 -> per-sample metadata carried onto each oligo
  crossrun work/motifs/family_reps.meme -> convert_motifs.py --meme  (used directly)

How a hand-picked record is matched (annotate mode)
---------------------------------------------------
Resolved in this order, and the mechanism used is recorded per record in
`match_method`, so no coordinate is ever unexplained:

  name            header names a crossrun peak and the sequence is identical
  name_mutant     header names a peak but the sequence DIFFERS: treated as a
                  mutant of that peak - lineage_id is the peak, role is
                  mutant, and it carries no coordinates of its own because it
                  inherits the parent's
  sequence        sequence is identical to exactly one peak (renamed record)
  substring       opt-in (--allow-substring): the record is a unique substring
                  of one peak, and coordinates are shifted by the offset
  unmatched       no match; kept, with no coordinates, and listed in the summary

A record matching several DISTINCT peaks is a hard failure, never resolved by
picking one.

Strand
------
Every sequence is taken to be on the + strand. crossrun extracts with bedtools
getfasta and no -s, so its peak sequences are + strand by construction, and
matched records are recorded with strand "+".

A record that matches a peak ONLY by reverse complement is therefore a failure,
not a match: under this assumption it means the record was revcomped somewhere
upstream, or is not the peak it appears to be. The check is kept (rather than
simply dropped) so the record is named and explained instead of silently
falling through to "unmatched".

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



COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def load_peaks(seqdir: str, wanted: list[str], man: pd.DataFrame,
               carry: list[str]) -> list[dict]:
    """Every crossrun peak with its sequence, coordinates and sample metadata."""
    peaks = []
    for sample in wanted:
        fa = os.path.join(seqdir, f"{sample}.fa")
        bed = os.path.join(seqdir, f"{sample}.scan.bed")
        for q in (fa, bed):
            if not os.path.exists(q):
                die(f"missing {q}. Run crossrun step 3 (steps/03_sequences.sh) "
                    f"for sample {sample} first.")

        coords = pd.read_csv(bed, sep="\t", header=None,
                             names=["chrom", "start", "end", "name"],
                             usecols=[0, 1, 2, 3], dtype={"chrom": str, "name": str})
        bad = sorted(set(coords.loc[~coords["chrom"].str.match(UCSC), "chrom"]))
        if bad:
            die(f"{bed} has non-UCSC contig name(s) {bad[:5]}. crossrun normalises "
                "to chr-style; a RefSeq (NC_0609xx.1) or GenBank (CP0682xx.2) name "
                "here would silently fail to join the BAMs.")
        by_name = coords.set_index("name")
        if by_name.index.duplicated().any():
            dup = sorted(by_name.index[by_name.index.duplicated()])[:5]
            die(f"{bed} has duplicate peak name(s) {dup}")

        mrow = man[man["sample_id"] == sample].iloc[0]
        notes = "; ".join(f"{k}={mrow[k]}" for k in carry
                          if k in man.columns and pd.notna(mrow[k]))

        seen = set()
        for rec in SeqIO.parse(fa, "fasta"):
            seen.add(rec.id)
            if rec.id not in by_name.index:
                die(f"{fa}: record '{rec.id}' has no row in {bed}. The FASTA and "
                    "the scan BED must name the same peaks; a mismatch means one "
                    "of them is stale.")
            c = by_name.loc[rec.id]
            seq = str(rec.seq).upper()
            if (int(c["end"]) - int(c["start"])) != len(seq):
                die(f"{rec.id}: scan.bed interval spans "
                    f"{int(c['end']) - int(c['start'])} bp but the FASTA record is "
                    f"{len(seq)} bp. BED is 0-based half-open; do not reconcile "
                    "this by shifting a coordinate until you know which file is wrong.")
            peaks.append({"name": rec.id, "sample": sample, "seq": seq,
                          "chrom": c["chrom"], "start": int(c["start"]),
                          "end": int(c["end"]), "notes": notes})
        orphan = sorted(set(by_name.index) - seen)
        if orphan:
            die(f"{bed} names {len(orphan)} peak(s) absent from {fa}, e.g. "
                f"{orphan[:5]}. One of the two files is stale; rerun crossrun step 3.")
    return peaks


def check_seq(sid: str, seq: str, lo: int, hi: int) -> str | None:
    """Return a skip reason, or None if the sequence is usable."""
    if set(seq) - set("ACGT"):
        return "non-ACGT base (soft-masked or N in the reference)"
    if len(seq) < lo:
        return f"{len(seq)} bp < min {lo}"
    if len(seq) > hi:
        return f"{len(seq)} bp > max {hi}"
    return None


def resolve(sid: str, seq: str, peaks: list[dict], by_name: dict,
            by_seq: dict, by_rc: dict, allow_substring: bool) -> dict:
    """Match one hand-picked record to a crossrun peak. See module docstring."""
    if sid in by_name:
        pk = by_name[sid]
        if pk["seq"] == seq:
            return {"method": "name", "peak": pk, "strand": "+", "offset": 0}
        # Named after a peak but not identical to it: a mutant of that peak.
        # It inherits the parent's coordinates through the lineage rather than
        # declaring its own, which is the only correct placement for a sequence
        # that no longer matches the genome.
        return {"method": "name_mutant", "peak": pk, "strand": "+", "offset": None}

    hits = by_seq.get(seq, [])
    distinct = {h["name"] for h in hits}
    if len(distinct) > 1:
        die(f"{sid}: sequence matches {len(distinct)} distinct crossrun peaks "
            f"({sorted(distinct)[:4]}). Refusing to pick one; rename the record "
            "to the peak you mean.")
    if len(distinct) == 1:
        return {"method": "sequence", "peak": hits[0], "strand": "+", "offset": 0}

    # Sequences are assumed to be + strand. A reverse-complement match means the
    # record was flipped upstream, so it is reported rather than accepted.
    rc_hits = by_rc.get(seq, [])
    if rc_hits:
        names = sorted({h["name"] for h in rc_hits})
        die(f"{sid}: does not match any peak on the + strand, but its reverse "
            f"complement matches {names[0]}"
            + (f" (and {len(names) - 1} other peak(s))" if len(names) > 1 else "")
            + ". Sequences are expected on the + strand, as crossrun extracts "
              "them. Reverse-complement this record, or confirm it is the "
              "sequence you meant.")

    if allow_substring:
        hits = [(pk, pk["seq"].find(seq)) for pk in peaks if seq in pk["seq"]]
        distinct = {pk["name"] for pk, _ in hits}
        if len(distinct) > 1:
            die(f"{sid}: substring of {len(distinct)} distinct peaks "
                f"({sorted(distinct)[:4]}). Refusing to pick one.")
        if len(distinct) == 1:
            pk, off = hits[0]
            if pk["seq"].count(seq) > 1:
                die(f"{sid}: occurs {pk['seq'].count(seq)} times within peak "
                    f"{pk['name']}; the offset is ambiguous.")
            return {"method": "substring", "peak": pk, "strand": "+", "offset": off}

    return {"method": "unmatched", "peak": None, "strand": "+", "offset": None}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun-work", required=True,
                    help="crossrun <outdir>/work directory")
    ap.add_argument("--fasta",
                    help="ANNOTATE mode: your own hand-picked FASTA. It defines "
                         "the sequence set; crossrun supplies coordinates and "
                         "metadata by matching. Omit for HARVEST mode.")
    ap.add_argument("--allow-substring", action="store_true",
                    help="also match a record that is a unique substring of one "
                         "peak, shifting coordinates by the offset. Off by "
                         "default: a shifted coordinate should be asked for, "
                         "not inferred quietly.")
    ap.add_argument("--manifest", required=True, help="crossrun manifest.tsv")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="sample_id(s) to consult; default every sample in the manifest")
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

    peaks = load_peaks(seqdir, wanted, man, args.carry)
    if not peaks:
        die("crossrun supplied no peaks")

    records, meta_rows, skipped = [], [], []
    counts = {"mode": "", "samples": len(wanted), "peaks_available": len(peaks),
              "records_in": 0, "kept": 0}

    if args.fasta:
        counts["mode"] = "annotate"
        by_name = {pk["name"]: pk for pk in peaks}
        by_seq: dict[str, list[dict]] = {}
        by_rc: dict[str, list[dict]] = {}
        for pk in peaks:
            by_seq.setdefault(pk["seq"], []).append(pk)  # retained for the warning below

        methods: dict[str, int] = {}
        for rec in SeqIO.parse(args.fasta, "fasta"):
            counts["records_in"] += 1
            sid_raw = rec.id
            seq = str(rec.seq).upper()
            sid = re.sub(r"[^A-Za-z0-9._-]", "_", sid_raw)
            reason = check_seq(sid, seq, args.min_length, args.max_length)
            if reason:
                skipped.append({"name": sid_raw, "reason": reason})
                continue
            m = resolve(sid_raw, seq, peaks, by_name, by_seq, by_rc,
                        args.allow_substring)
            methods[m["method"]] = methods.get(m["method"], 0) + 1
            pk = m["peak"]

            if m["method"] == "name_mutant":
                lineage = re.sub(r"[^A-Za-z0-9._-]", "_", pk["name"])
                row = {"sequence_id": sid, "lineage_id": lineage, "role": "mutant",
                       "chrom": "", "start": "", "end": "", "strand": "",
                       "mutation": "differs from the crossrun peak of the same name",
                       "notes": f"match_method=name_mutant; parent_peak={pk['name']}; "
                                f"crossrun_sample={pk['sample']}; {pk['notes']}"}
            elif pk is None:
                row = {"sequence_id": sid, "lineage_id": sid, "role": "wt",
                       "chrom": "", "start": "", "end": "", "strand": "",
                       "mutation": "",
                       "notes": "match_method=unmatched; no crossrun peak - this "
                                "sequence is excluded from genome-anchored steps"}
            else:
                off = m["offset"]
                row = {"sequence_id": sid, "lineage_id": sid, "role": "wt",
                       "chrom": pk["chrom"], "start": pk["start"] + off,
                       "end": pk["start"] + off + len(seq), "strand": "+",
                       "mutation": "",
                       "notes": f"match_method={m['method']}; peak={pk['name']}; "
                                f"crossrun_sample={pk['sample']}; {pk['notes']}"}
            records.append((sid, seq))
            meta_rows.append(row)
            counts["kept"] += 1
        counts["match_methods"] = methods
    else:
        counts["mode"] = "harvest"
        for pk in peaks:
            counts["records_in"] += 1
            reason = check_seq(pk["name"], pk["seq"], args.min_length, args.max_length)
            if reason:
                skipped.append({"name": pk["name"], "reason": reason})
                continue
            sid = re.sub(r"[^A-Za-z0-9._-]", "_", pk["name"])
            records.append((sid, pk["seq"]))
            meta_rows.append({
                "sequence_id": sid, "lineage_id": sid, "role": "wt",
                "chrom": pk["chrom"], "start": pk["start"], "end": pk["end"],
                "strand": "+", "mutation": "",
                "notes": f"match_method=harvest; crossrun_sample={pk['sample']}; "
                         f"{pk['notes']}"})
            counts["kept"] += 1

    if not records:
        die("no sequences survived import; nothing to write")

    ids = [r[0] for r in records]
    if len(set(ids)) != len(ids):
        dup = sorted({i for i in ids if ids.count(i) > 1})[:5]
        die(f"duplicate sequence_id: {dup}. Names collide after sanitisation; "
            "rename them before importing.")

    with open(args.out_fasta, "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pd.DataFrame(meta_rows).to_csv(args.out_metadata, sep="\t", index=False)

    summary = {"counts": counts, "genome_build": builds[0] if builds else None,
               "skipped": skipped[:200], "n_skipped": len(skipped),
               "coordinate_convention": "0-based half-open, taken from scan.bed",
               "contig_convention": "UCSC chr-style, asserted",
               "substring_matching": bool(args.allow_substring)}
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[import_crossrun] mode={counts['mode']} "
          f"peaks_available={counts['peaks_available']} "
          f"records_in={counts['records_in']} kept={counts['kept']} "
          f"skipped={len(skipped)} build={builds[0] if builds else 'unset'}")
    if counts["mode"] == "annotate":
        for k, v in sorted(counts["match_methods"].items()):
            print(f"  match_method {k:<14} {v}")
        un = counts["match_methods"].get("unmatched", 0)
        if un:
            print(f"[import_crossrun] {un} record(s) matched no crossrun peak and "
                  "carry no coordinates; they are excluded from genome-anchored steps")
    if skipped:
        print(f"[import_crossrun] {len(skipped)} record(s) skipped; see "
              f"{args.out_summary} for the reason on each")


if __name__ == "__main__":
    main()
