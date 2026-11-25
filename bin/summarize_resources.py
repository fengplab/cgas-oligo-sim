#!/usr/bin/env python3
"""Turn a Nextflow execution trace into per-process resource recommendations.

Reads the raw trace written by the profiling block in nextflow.config and
reports, per process, what was REQUESTED against what was actually USED, then
emits a paste-ready Groovy block.

Two deliberate choices, because getting them wrong is expensive in opposite
directions:

  * Memory is recommended from the MAXIMUM observed peak_rss, not a quantile.
    A quantile saves allocation but OOM-kills the tail, and an OOM-killed task
    at hour 20 of a 51,000-task campaign costs more than the memory ever saved.
  * Only COMPLETED tasks inform the recommendation. A failed task's peak_rss
    records where it was when it died, not what the work needs, and it can
    mislead in either direction: an OOM kill pegs peak_rss at the limit and
    inflates the next request, while an early unrelated failure truncates it
    low and deflates it. Neither is a measurement of the requirement.

Failed and retried tasks are reported separately rather than dropped silently.
"""

import argparse
import glob
import math
import sys
from collections import defaultdict

# Raw-mode trace units, verified against Nextflow 25.10.2:
#   peak_rss / peak_vmem / rchar / wchar : bytes
#   duration / realtime                  : milliseconds
#   %cpu                                 : bare number, may exceed 100
#   cpus                                 : integer; memory: bytes, or "-" if unset
MISSING = {"-", "", "null", None}


def die(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(1)


def as_float(v):
    if v in MISSING:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def human_bytes(n) -> str:
    if n is None:
        return "-"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
        n /= 1024.0
    return f"{n:.1f} TB"


def read_trace(paths):
    rows, total = [], 0
    for path in paths:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            needed = {"process", "status", "peak_rss", "%cpu", "realtime", "cpus"}
            missing = needed - set(header)
            if missing:
                die(f"{path} lacks trace fields {sorted(missing)}. It was probably "
                    "written before the profiling block was added, or with a "
                    "different `fields` list.")
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                total += 1
                rows.append(dict(zip(header, line.split("\t"))))
    if not rows:
        die("no task rows in the supplied trace file(s)")
    return rows, total


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", nargs="+", required=True,
                    help="execution_trace_*.txt; globs accepted")
    ap.add_argument("--memory-headroom", type=float, default=1.5,
                    help="multiply max observed peak_rss by this for the "
                         "recommendation (default 1.5)")
    ap.add_argument("--min-memory-gb", type=float, default=1.0,
                    help="never recommend less than this (default 1)")
    ap.add_argument("--cpu-floor", type=float, default=0.6,
                    help="flag a process whose mean %%cpu per requested core "
                         "falls below this fraction (default 0.6)")
    args = ap.parse_args()

    paths = []
    for pattern in args.trace:
        hits = sorted(glob.glob(pattern))
        paths.extend(hits if hits else [pattern])
    rows, n_rows = read_trace(paths)

    by_proc = defaultdict(list)
    for r in rows:
        by_proc[r["process"]].append(r)

    completed = [r for r in rows if r.get("status") == "COMPLETED"]
    other = [r for r in rows if r.get("status") != "COMPLETED"]
    retried = [r for r in rows if as_float(r.get("attempt")) and
               as_float(r["attempt"]) > 1]

    print(f"trace files        : {len(paths)}")
    print(f"task rows read     : {n_rows}")
    print(f"  COMPLETED        : {len(completed)}")
    print(f"  other status     : {len(other)}")
    print(f"  attempt > 1      : {len(retried)}")
    if other:
        seen = defaultdict(int)
        for r in other:
            seen[(r["process"], r.get("status"), r.get("exit"))] += 1
        print("\nnot COMPLETED (excluded from recommendations):")
        for (proc, st, ex), n in sorted(seen.items()):
            flag = "  <- exit 137/140 often means the scheduler killed it (OOM/limits)" \
                   if ex in {"137", "140"} else ""
            print(f"  {proc:<22} {st:<12} exit={ex:<5} x{n}{flag}")

    print("\n" + "=" * 94)
    print(f"{'process':<22}{'n':>5}{'req cpus':>9}{'mean %cpu':>11}"
          f"{'req mem':>11}{'max RSS':>11}{'use':>7}{'max real':>11}")
    print("=" * 94)

    recs = {}
    for proc in sorted(by_proc):
        ok = [r for r in by_proc[proc] if r.get("status") == "COMPLETED"]
        if not ok:
            print(f"{proc:<22}{len(by_proc[proc]):>5}   no COMPLETED task; cannot recommend")
            continue
        rss = [v for v in (as_float(r.get("peak_rss")) for r in ok) if v is not None]
        cpu = [v for v in (as_float(r.get("%cpu")) for r in ok) if v is not None]
        real = [v for v in (as_float(r.get("realtime")) for r in ok) if v is not None]
        req_cpus = next((as_float(r.get("cpus")) for r in ok
                         if as_float(r.get("cpus"))), None)
        req_mem = next((as_float(r.get("memory")) for r in ok
                        if as_float(r.get("memory"))), None)

        max_rss = max(rss) if rss else None
        mean_cpu = sum(cpu) / len(cpu) if cpu else None
        max_real = max(real) if real else None
        frac = f"{max_rss / req_mem:.0%}" if (max_rss and req_mem) else "-"

        print(f"{proc:<22}{len(ok):>5}"
              f"{(f'{req_cpus:.0f}' if req_cpus else '-'):>9}"
              f"{(f'{mean_cpu:.0f}%' if mean_cpu is not None else '-'):>11}"
              f"{human_bytes(req_mem):>11}{human_bytes(max_rss):>11}{frac:>7}"
              f"{(f'{max_real/1000:.1f}s' if max_real is not None else '-'):>11}")

        if max_rss is not None:
            gb = max(args.min_memory_gb,
                     math.ceil(max_rss * args.memory_headroom / 2**30 * 2) / 2)
            cpus_rec = None
            if mean_cpu is not None and req_cpus:
                used = mean_cpu / 100.0
                if used < args.cpu_floor * req_cpus:
                    candidate = max(1, math.ceil(used))
                    # Only surface it if it actually differs from the request;
                    # "reduce 1 cpu to 1 cpu" is noise, and noise in a
                    # recommendation block is how real findings get skipped.
                    if candidate < req_cpus:
                        cpus_rec = candidate
            recs[proc] = (gb, cpus_rec, mean_cpu, req_cpus)

    print("=" * 94)
    print("'use' is max RSS over requested memory. Under ~25% is over-allocated;")
    print("over 100% means the task survived only because the scheduler was lenient.")

    if not recs:
        print("\nNo process had a COMPLETED task with a peak_rss, so there is nothing "
              "to recommend.")
        return

    print("\n" + "-" * 94)
    print(f"Recommended, from max observed RSS x {args.memory_headroom}:")
    print("-" * 94)
    print("process {")
    for proc, (gb, cpus_rec, mean_cpu, req_cpus) in sorted(recs.items()):
        print(f"    withName: '{proc}' {{")
        print(f"        memory = '{gb:g} GB'")
        if cpus_rec is not None:
            print(f"        cpus   = {cpus_rec}"
                  f"        // mean {mean_cpu:.0f}% of {req_cpus:.0f} requested")
        print("    }")
    print("}")
    print("\nCheck before pasting: a process whose label sets clusterOptions must")
    print("keep --ntasks and --cpus-per-task, and no `time` directive may appear.")


if __name__ == "__main__":
    main()
