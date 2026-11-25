#!/usr/bin/env python3
"""Combine the per-arm site calls into one table plus a cross-arm comparison.

Arms are called separately, because call_sites.py pools over replicate labels
and would otherwise average across different force fields and structures. This
script puts them back together WITH the arm as a column, and reports where they
agree and where they do not.

What it deliberately does not do: no significance test between arms. The arms
differ in force field, reference structure or an imposed motif term, so a
difference in occupancy is a difference in model, not an measurement of the same
quantity under two conditions. Treating a cross-arm gap as an effect size with a
p-value would be reading a modelling choice as a result. The tables report the
difference and leave the interpretation to the methods.
"""

import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")           # set before pyplot: there is no display here
import matplotlib.pyplot as plt
import pandas as pd

SITE_KEY = ["sequence_id", "site_id"]


def die(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(1)


def read_labelled(paths, what):
    """Read arm=path pairs into one frame with an `arm` column."""
    frames = []
    for spec in paths:
        if "=" not in spec:
            die(f"--{what} wants arm=path, got {spec!r}")
        arm, path = spec.split("=", 1)
        if not arm:
            die(f"--{what} has an empty arm label: {spec!r}")
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            print(f"warning: {path} (arm {arm}) has no rows", file=sys.stderr)
        df.insert(0, "arm", arm)
        frames.append(df)
    if not frames:
        die(f"no --{what} given")
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", nargs="+", required=True,
                    help="arm=path to that arm's sites.tsv")
    ap.add_argument("--profiles", nargs="+", required=True,
                    help="arm=path to that arm's occupancy profile")
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--out-comparison", required=True)
    ap.add_argument("--out-figures", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    sites = read_labelled(args.sites, "sites")
    profiles = read_labelled(args.profiles, "profiles")

    for col in SITE_KEY:
        if col not in sites.columns:
            die(f"combined sites table lacks {col}")
    sites.to_csv(args.out_sites, sep="\t", index=False)

    arms = sorted(sites["arm"].unique())
    n_arms = len(arms)

    # Which arms called each site. A site called in one arm only is the
    # interesting case: it says the call depends on the model.
    presence = (sites.assign(called=1)
                .pivot_table(index=SITE_KEY, columns="arm", values="called",
                             aggfunc="max", fill_value=0))
    for a in arms:
        if a not in presence.columns:
            presence[a] = 0
    presence = presence[arms]
    presence["n_arms_calling"] = presence.sum(axis=1)
    presence["arms_calling"] = presence[arms].apply(
        lambda r: ",".join([a for a in arms if r[a] == 1]), axis=1)

    occ_col = "occupancy" if "occupancy" in sites.columns else None
    if occ_col:
        wide = sites.pivot_table(index=SITE_KEY, columns="arm", values=occ_col,
                                 aggfunc="max")
        wide.columns = [f"occupancy_{c}" for c in wide.columns]
        presence = presence.join(wide)
        if n_arms >= 2:
            cols = [f"occupancy_{a}" for a in arms if f"occupancy_{a}" in presence]
            if len(cols) >= 2:
                # max-min skips NaN, so a site called in ONE arm would get a
                # range of 0.0 - which reads as perfect agreement between arms
                # when in fact there is nothing to compare. Only define it where
                # at least two arms actually reported an occupancy.
                have = presence[cols].notna().sum(axis=1)
                rng = presence[cols].max(axis=1) - presence[cols].min(axis=1)
                presence["occupancy_range"] = rng.where(have >= 2)

    presence.reset_index().to_csv(args.out_comparison, sep="\t", index=False)

    shared = int((presence["n_arms_calling"] == n_arms).sum())
    unique = int((presence["n_arms_calling"] == 1).sum())

    # ---- figures ----
    pos_col = "position" if "position" in profiles.columns else None
    occ_prof = "occupancy" if "occupancy" in profiles.columns else None

    n_panels = 2 if (pos_col and occ_prof) else 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 4 * n_panels),
                             squeeze=False)

    ax = axes[0][0]
    counts = sites.groupby("arm").size().reindex(arms, fill_value=0)
    ax.bar(range(len(arms)), counts.to_numpy())
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(arms, rotation=20, ha="right")
    ax.set_ylabel("sites called")
    ax.set_title(f"Sites per arm  (shared by all {n_arms}: {shared}; "
                 f"called in one arm only: {unique})")

    if n_panels == 2:
        ax = axes[1][0]
        for a in arms:
            sub = (profiles[profiles["arm"] == a]
                   .groupby(pos_col, as_index=False)[occ_prof].mean()
                   .sort_values(pos_col))
            ax.plot(sub[pos_col], sub[occ_prof], label=a, linewidth=1.2)
        ax.set_xlabel("position (bp)")
        ax.set_ylabel("mean occupancy")
        ax.set_title("Occupancy profile by arm")
        ax.legend(fontsize="small")

    fig.tight_layout()
    fig.savefig(args.out_figures, dpi=150)
    plt.close(fig)

    with open(args.out_summary, "w") as fh:
        json.dump({"n_arms": n_arms, "arms": arms,
                   "n_sites_total_rows": int(len(sites)),
                   "n_distinct_sites": int(len(presence)),
                   "sites_called_in_all_arms": shared,
                   "sites_called_in_one_arm_only": unique,
                   "note": "no cross-arm significance test; arms differ in "
                           "model, not in condition"}, fh, indent=2)

    print(f"[cross_arm] arms={n_arms} distinct sites={len(presence)} "
          f"shared={shared} single-arm={unique}")


if __name__ == "__main__":
    main()
