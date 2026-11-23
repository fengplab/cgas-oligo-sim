#!/usr/bin/env python3
"""Estimate the compute budget for a run configuration.

The trajectory COUNT is exact arithmetic. The per-trajectory COST is not
known until you measure it: no GENESIS benchmark on this cluster exists, and
a throughput figure taken from anywhere else would be a guess dressed as a
number. Measure one trajectory, pass --min-per-traj, and the rest is
multiplication.

Arms follow the baseline + one-factor-at-a-time design: a baseline arm, one
arm per extra partner, one arm per motif tested with PWMcos, one arm per
extra structure, plus a multi-basin arm when two or more structures are
given. Mutants are NOT an arm: they are sequences in the panel.
"""
from __future__ import annotations
import argparse


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--oligos", type=int, default=50, help="sequences in the panel, parents + mutants")
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--shuffles", type=int, default=20, help="null shuffles per sequence")
    ap.add_argument("--null-replicates", type=int, default=12)
    ap.add_argument("--partners", type=int, default=1, help="extra partner arms beyond baseline")
    ap.add_argument("--motifs", type=int, default=10, help="PWMcos arms")
    ap.add_argument("--structures", type=int, default=1)
    ap.add_argument("--min-per-traj", type=float, default=None,
                    help="MEASURED wall-minutes for one trajectory")
    ap.add_argument("--concurrent", type=int, default=100, help="trajectories running at once")
    args = ap.parse_args()

    arms = 1 + args.partners + args.motifs + max(0, args.structures - 1)
    if args.structures >= 2:
        arms += 1
    obs = args.oligos * args.replicates
    null = args.oligos * args.shuffles * args.null_replicates
    per_arm = obs + null
    total = arms * per_arm

    print(f"arms                     {arms}")
    print(f"  baseline               1")
    print(f"  partner arms           {args.partners}")
    print(f"  PWMcos arms            {args.motifs}")
    print(f"  extra structure arms   {max(0, args.structures - 1)}")
    print(f"  multi-basin arm        {1 if args.structures >= 2 else 0}")
    print()
    print(f"per arm, observed        {obs:>9,}  ({args.oligos} x {args.replicates})")
    print(f"per arm, null            {null:>9,}  ({args.oligos} x {args.shuffles} x {args.null_replicates})")
    print(f"per arm, total           {per_arm:>9,}   null is {100*null/per_arm:.0f}% of it")
    print()
    print(f"TOTAL TRAJECTORIES       {total:>9,}")

    if args.min_per_traj is None:
        print()
        print("No --min-per-traj given, so no wall-time estimate.")
        print("Measure it:  time mpirun -np 4 atdyn md.inp   on one representative system.")
        for m in (1, 5, 15, 60):
            cpu_h = total * m / 60
            wall_d = cpu_h / args.concurrent / 24
            print(f"  if {m:>3} min/traj -> {cpu_h:>12,.0f} core-hours, "
                  f"{wall_d:>7.1f} days wall at {args.concurrent} concurrent")
    else:
        cpu_h = total * args.min_per_traj / 60
        wall_d = cpu_h / args.concurrent / 24
        print(f"at {args.min_per_traj} min/traj -> {cpu_h:,.0f} core-hours, "
              f"{wall_d:.1f} days wall at {args.concurrent} concurrent")


if __name__ == "__main__":
    main()
