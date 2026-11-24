#!/usr/bin/env python3
"""Assert the stored span equals the sequence length, whatever the input said.

That invariant is what the rest of the pipeline relies on: every module
downstream treats metadata coordinates as 0-based half-open, so if conversion
is wrong here it is wrong everywhere, silently.
"""
import json
import sys

import pandas as pd

conv = sys.argv[1]
md = pd.read_csv(f"tests/out/coords_{conv}.tsv", sep="\t")

fa, name = {}, None
for line in open(f"tests/out/coords_{conv}.built.fa"):
    if line.startswith(">"):
        name = line[1:].strip()
    else:
        fa[name] = fa.get(name, "") + line.strip()

for _, r in md.iterrows():
    span = int(r.end) - int(r.start)
    n = len(fa[r.sequence_id])
    if span != n:
        sys.exit(f"FAIL: {r.sequence_id} stored span {span} != length {n}")

ev = json.load(open(f"tests/out/coords_{conv}.json"))["conventions"]
if ev["header_coords_mode"] != "auto":
    sys.exit("FAIL: summary should record that detection ran")
if not ev["header_coords_evidence"]["decision"].startswith("unanimous"):
    sys.exit(f"FAIL: expected a unanimous decision, got "
             f"{ev['header_coords_evidence']['decision']!r}")
print(f"  OK {conv}: detected, stored spans match sequence lengths")
