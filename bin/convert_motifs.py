#!/usr/bin/env python3
"""Convert MEME PWMs or IUPAC consensus strings into PFMs for GENESIS PWMcos.

Accepts
-------
  --meme FILE    MEME Suite output or minimal-MEME format (Bio.motifs)
  --iupac FILE   TSV: motif_id, consensus   (IUPAC degenerate DNA)

Emits one PFM per motif plus a manifest.

PWMcos is an ENERGY TERM, not an annotation
-------------------------------------------
A PFM handed to GENESIS via `aa_2_cg.jl --pwmcos -p <pfm>` becomes part of the
force field: it imposes the sequence preference it encodes. Runs using it
cannot answer "does cGAS prefer this motif" - that is circular. They answer
"what would this motif predict, if cGAS read it". Keep PWMcos-on and
PWMcos-off runs in separate, labelled arms and never pool them.

# UNVERIFIED: the exact GENESIS PFM layout. The published tutorial edits a
# JASPAR matrix by hand to add CHAIN_A / CHAIN_B index lines and to replace
# zeros with pseudocounts, but does not state a formal grammar. What is written
# here is JASPAR-style counts plus those two chain lines. Check one output
# against `aa_2_cg.jl --pwmcos` on a known case before trusting a production run.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
BASES = "ACGT"


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def iupac_to_counts(consensus: str, depth: int = 100) -> np.ndarray:
    """Expand a degenerate consensus into a count matrix (4 x L)."""
    cons = consensus.strip().upper()
    mat = np.zeros((4, len(cons)))
    for j, ch in enumerate(cons):
        if ch not in IUPAC:
            die(f"'{ch}' is not an IUPAC DNA code (in '{cons}')")
        allowed = IUPAC[ch]
        for b in allowed:
            mat[BASES.index(b), j] = depth / len(allowed)
    return mat


def read_meme(path: str) -> list[tuple[str, np.ndarray]]:
    from Bio import motifs
    last = None
    for fmt in ("meme", "minimal"):
        try:
            with open(path) as fh:
                recs = motifs.parse(fh, fmt)
            out = []
            for m in recs:
                counts = np.array([list(m.counts[b]) for b in BASES], dtype=float)
                name = m.name or m.consensus
                out.append((str(name), counts))
            if out:
                return out
        except Exception as exc:          # try the next dialect
            last = exc
    die(f"could not parse {path} as MEME or minimal-MEME format: {last}")


def write_pfm(path: str, motif_id: str, counts: np.ndarray,
              pseudocount: float, chain_a: str, chain_b: str) -> dict:
    """JASPAR-style counts + GENESIS chain index lines. See UNVERIFIED note."""
    n_zero = int((counts == 0).sum())
    counts = counts + pseudocount          # zeros break the log-odds in PWMcos
    with open(path, "w") as fh:
        fh.write(f"> {motif_id}\n")
        fh.write(f"CHAIN_A {chain_a}\n")
        fh.write(f"CHAIN_B {chain_b}\n")
        for i, b in enumerate(BASES):
            vals = " ".join(f"{v:10.4f}" for v in counts[i])
            fh.write(f"{b} [ {vals} ]\n")
    return {"motif_id": motif_id, "length": int(counts.shape[1]),
            "zeros_replaced": n_zero, "pseudocount": pseudocount,
            "pfm_path": os.path.abspath(path)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--meme")
    src.add_argument("--iupac")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pseudocount", type=float, default=0.25,
                    help="added to every cell; zeros make the PWMcos log-odds "
                         "undefined, so this is not cosmetic")
    ap.add_argument("--chain-a", default="1",
                    help="DNA chain index of the sense strand, as GENESIS numbers "
                         "chains in the .top [molecules] block")
    ap.add_argument("--chain-b", default="2", help="antisense strand chain index")
    ap.add_argument("--out-manifest", required=True)
    args = ap.parse_args()

    if args.pseudocount <= 0:
        die(f"--pseudocount must be > 0, got {args.pseudocount}")
    os.makedirs(args.outdir, exist_ok=True)

    if args.meme:
        parsed = read_meme(args.meme)
        source = os.path.abspath(args.meme)
        kind = "meme"
    else:
        parsed, source, kind = [], os.path.abspath(args.iupac), "iupac"
        with open(args.iupac) as fh:
            for ln, line in enumerate(fh, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    die(f"{args.iupac}:{ln} expected 'motif_id<TAB>consensus'")
                if parts[0].lower() in ("motif_id", "id", "name"):
                    continue
                parsed.append((parts[0], iupac_to_counts(parts[1])))

    if not parsed:
        die(f"no motifs parsed from {source}")

    seen, manifest = set(), []
    for motif_id, counts in parsed:
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in motif_id)
        if safe in seen:
            die(f"duplicate motif id after sanitisation: {safe}")
        seen.add(safe)
        if counts.shape[0] != 4:
            die(f"motif {motif_id}: expected 4 rows (ACGT), got {counts.shape[0]}")
        if counts.shape[1] < 4:
            die(f"motif {motif_id}: length {counts.shape[1]} is too short to be a "
                "usable recognition motif")
        if not np.isfinite(counts).all() or (counts < 0).any():
            die(f"motif {motif_id}: counts contain negative or non-finite values")
        manifest.append(write_pfm(os.path.join(args.outdir, f"{safe}.pfm"),
                                  motif_id, counts, args.pseudocount,
                                  args.chain_a, args.chain_b))

    with open(args.out_manifest, "w") as fh:
        json.dump({"source": source, "format": kind,
                   "pseudocount": args.pseudocount,
                   "chain_a": args.chain_a, "chain_b": args.chain_b,
                   "n_motifs": len(manifest), "motifs": manifest,
                   "warning": "PWMcos PFMs are force-field terms. Runs using them "
                              "impose the motif and cannot test for it."}, fh, indent=2)

    print(f"[convert_motifs] {kind}: {len(manifest)} motif(s) -> {args.outdir}")
    for m in manifest:
        print(f"  {m['motif_id']:<24} len={m['length']:<4} "
              f"zeros_replaced={m['zeros_replaced']}")


if __name__ == "__main__":
    main()
