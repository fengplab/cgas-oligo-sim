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


def budget(n_traj, min_per_traj, omp_threads, concurrent):
    """Core-hours and wall days for n_traj trajectories.

    These scale differently and must not be derived from each other. Core-hours
    count the threads each trajectory occupies, because a serial atdyn with
    OMP_NUM_THREADS=n holds n cores for its whole wall time. Wall time does NOT:
    `concurrent` counts trajectories in flight, not cores, so folding the thread
    count into it would inflate the schedule by exactly that factor.
    """
    traj_minutes = n_traj * min_per_traj
    cpu_h = traj_minutes * omp_threads / 60.0
    wall_d = traj_minutes / concurrent / 60.0 / 24.0
    return cpu_h, wall_d


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
    ap.add_argument("--omp-threads", type=int, default=2,
                    help="OpenMP threads per trajectory; must match params.cg_omp_threads. "
                         "atdyn is serial + OpenMP, so this is the core count per trajectory")
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
        print("Measure it:  OMP_NUM_THREADS=%d time atdyn md.inp   on one representative"
              % args.omp_threads)
        print("system. atdyn is built --disable-mpi, so it is one serial process with")
        print("OpenMP threads; there are no ranks.")
        for m in (1, 5, 15, 60):
            cpu_h, wall_d = budget(total, m, args.omp_threads, args.concurrent)
            print(f"  if {m:>3} min/traj -> {cpu_h:>12,.0f} core-hours, "
                  f"{wall_d:>7.1f} days wall at {args.concurrent} concurrent")
    else:
        cpu_h, wall_d = budget(total, args.min_per_traj, args.omp_threads,
                               args.concurrent)
        print(f"at {args.min_per_traj} min/traj -> {cpu_h:,.0f} core-hours, "
              f"{wall_d:.1f} days wall at {args.concurrent} concurrent "
              f"({args.omp_threads} thread(s) each)")


if __name__ == "__main__":
    main()
