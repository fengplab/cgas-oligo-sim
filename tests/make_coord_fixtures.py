#!/usr/bin/env python3
"""FASTA fixtures for coordinate-convention detection.

Exactly one convention can fit a record, because end-start and end-start+1
cannot both equal the same length. So detection is a unanimity test, and these
four files cover its outcomes: unanimous 1-closed, unanimous 0-half-open, a
file that mixes the two, and a file whose span matches neither reading.
"""
import random

random.seed(11)
OUT = "tests/out"


def seq(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def write(path: str, rows: list[tuple[int, int]], closed: bool) -> None:
    with open(path, "w") as fh:
        for i, (lo, n) in enumerate(rows, 1):
            end = lo + n - 1 if closed else lo + n
            fh.write(f">rank{i} chr1:{lo}-{end} summit={lo + n // 2}\n")
            fh.write(seq(n) + "\n")


write(f"{OUT}/coords_1closed.fa", [(1000, 60), (5000, 80)], closed=True)
write(f"{OUT}/coords_0half.fa", [(1000, 60), (5000, 80)], closed=False)

# One record of each convention: the file mixes them.
with open(f"{OUT}/coords_mixed.fa", "w") as fh:
    fh.write(">a chr1:100-160 summit=130\n" + "A" * 60 + "\n")
    fh.write(">b chr1:200-260 summit=230\n" + "C" * 61 + "\n")

# Span matches neither reading: a real header/sequence disagreement.
with open(f"{OUT}/coords_neither.fa", "w") as fh:
    fh.write(">a chr1:100-160 summit=130\n" + "A" * 75 + "\n")

print("  wrote 4 coordinate fixtures")
