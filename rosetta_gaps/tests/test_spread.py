#!/usr/bin/env python3
"""Unit-test loop_spread against answers that can be worked out by hand.

The end-to-end fixture can only produce models whose loops differ slightly:
a larger random displacement of a loop's C-alphas breaks its peptide bonds, and
check_models.py then rightly rejects the model. Rather than loosen that bond
check to make a fixture pass, the spread arithmetic is tested directly on
constructed coordinates where the expected value is exact.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import gemmi

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "check_models", HERE.parent / "bin" / "check_models.py")
cm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cm)

FAILURES: list[str] = []


def check(name: str, got: float, want: float, tol: float = 1e-6) -> None:
    ok = got is not None and abs(got - want) < tol
    print(f"{'ok  ' if ok else 'FAIL'} {name}: got {got!r}, want {want!r}")
    if not ok:
        FAILURES.append(name)


def model(positions: list[tuple[float, float, float]]):
    """A (path, structure, ca-list) triple shaped like check_models expects."""
    return (Path("synthetic"), None,
            [gemmi.Position(*p) for p in positions])


# Six anchor points off a line, so superposition is well conditioned, plus
# three loop points.
FIXED = [(0, 0, 0), (3.8, 0.4, 0), (7.6, 1.6, 0.5),
         (11.4, 3.6, 1.5), (15.2, 6.4, 3.0), (19.0, 10.0, 5.0)]
LOOP = [(5.0, -4.0, 0.0), (8.0, -5.0, 1.0), (11.0, -4.0, 2.0)]
fixed_idx = list(range(len(FIXED)))
loop_idx = [len(FIXED) + i for i in range(len(LOOP))]

a = model(FIXED + LOOP)

# 1. A model identical to the first: spread is exactly zero, and the identical
#    fixed sets must not make the superposition return NaN.
d = cm.loop_spread([a, model(FIXED + LOOP)], fixed_idx, loop_idx)
check("identical models", d[0] if d else None, 0.0)

# 2. Every loop atom moved 2.0 A in one direction: RMSD over the loop is 2.0.
moved = [(x, y + 2.0, z) for x, y, z in LOOP]
d = cm.loop_spread([a, model(FIXED + moved)], fixed_idx, loop_idx)
check("loop displaced 2.0 A", d[0] if d else None, 2.0)

# 3. One loop atom of three moved 3.0 A: RMSD = sqrt((9+0+0)/3) = sqrt(3).
one = [(LOOP[0][0], LOOP[0][1] + 3.0, LOOP[0][2]), LOOP[1], LOOP[2]]
d = cm.loop_spread([a, model(FIXED + one)], fixed_idx, loop_idx)
check("one of three moved 3.0 A", d[0] if d else None, math.sqrt(3.0), 1e-6)

# 4. The WHOLE model rigidly translated: spread must be zero, because the
#    comparison superposes on the fixed residues first. If it were not, every
#    pair of Rosetta outputs would show a spread that is only a difference of
#    coordinate frame.
shift = [(x + 17.0, y - 9.0, z + 4.5) for x, y, z in FIXED + LOOP]
d = cm.loop_spread([a, model(shift)], fixed_idx, loop_idx)
check("whole model translated", d[0] if d else None, 0.0, 1e-4)

# 5. The whole model rigidly ROTATED: same requirement.
rot = gemmi.Transform()
ang = 0.7
rot.mat.fromlist([[math.cos(ang), -math.sin(ang), 0.0],
                  [math.sin(ang), math.cos(ang), 0.0],
                  [0.0, 0.0, 1.0]])
turned = [tuple(rot.apply(gemmi.Position(*p)).tolist()) for p in FIXED + LOOP]
d = cm.loop_spread([a, model(turned)], fixed_idx, loop_idx)
check("whole model rotated", d[0] if d else None, 0.0, 1e-4)

# 6. Three models give three pairs.
d = cm.loop_spread([a, model(FIXED + moved), model(FIXED + one)],
                   fixed_idx, loop_idx)
check("three models -> three pairs", float(len(d)), 3.0)

if FAILURES:
    sys.exit(f"\n{len(FAILURES)} spread test(s) failed: {FAILURES}")
print("\nall spread tests passed")
