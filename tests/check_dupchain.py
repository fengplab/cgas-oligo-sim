#!/usr/bin/env python3
"""The protein must survive stripping when chain names repeat."""
import sys

import gemmi
sys.path.insert(0, "tests")
from make_dupchain_fixture import build, ONE3  # noqa: E402

AA = set(ONE3.values())
st = build()
st.setup_entities()
st.remove_ligands_and_waters()

# the fix, as prepare_structures.py applies it
st.remove_empty_chains()
st.merge_chain_parts()

names = [c.name for c in st[0]]
if len(names) != len(set(names)):
    sys.exit(f"FAIL: chain names still repeat after cleanup: {names}")

drop_idx = [i for i, c in enumerate(st[0])
            if not ({r.name for r in c} & AA)]
for i in reversed(drop_idx):
    del st[0][i]

left = [(c.name, len(c)) for c in st[0]]
if not any({r.name for r in c} & AA for c in st[0]):
    sys.exit(f"FAIL: the protein chain was deleted; remaining {left}")
if [n for n, _ in left] != ["A"]:
    sys.exit(f"FAIL: expected only chain A to remain, got {left}")
print(f"  OK duplicate chain names handled: protein survived, remaining {left}")
