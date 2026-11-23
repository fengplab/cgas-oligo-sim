#!/usr/bin/env python3
"""Static figures for called occupancy sites.

Every figure ships a CSV companion (the table-view twin), so no value is
reachable only through the image. Colour scales are fixed across panels: a
per-panel autoscale makes noise look like signal when panels are compared.
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator

# --- palette (light mode) -------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
S1 = "#2a78d6"          # categorical slot 1 - observed
SEQ_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIV_LOW, DIV_MID, DIV_HIGH = "#2a78d6", "#f0efec", "#d03b3b"

BLUES = LinearSegmentedColormap.from_list("seq_blue", SEQ_RAMP)
DIVERGE = LinearSegmentedColormap.from_list("div_br", [DIV_LOW, DIV_MID, DIV_HIGH])

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK2,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 1.0,
    "grid.linestyle": "-",          # never dashed
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.facecolor": SURFACE,
})


def tidy(ax) -> None:
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


def fig_occupancy_track(prof: pd.DataFrame, sites: pd.DataFrame, seq: str,
                        outdir: str, ymax: float) -> None:
    """Figure 1 - occupancy along the sequence, against the shuffle null."""
    p = prof[prof["sequence_id"] == seq].sort_values("position")
    s = sites[sites["sequence_id"] == seq] if len(sites) else sites
    thr = float(p["null_threshold"].iloc[0])

    fig, ax = plt.subplots(figsize=(10, 3.4))
    # Null region: everything at or below the threshold is indistinguishable
    # from shuffled sequence. Shown as a region, not a series.
    ax.axhspan(0, thr, color=GRID, alpha=0.55, lw=0, zorder=0)
    ax.axhline(thr, color=MUTED, lw=1.0, zorder=1,
               label=f"shuffle null p{p['null_percentile'].iloc[0]:g} = {thr:.3f}"
               if "null_percentile" in p.columns else f"shuffle null = {thr:.3f}")
    ax.plot(p["position"], p["occupancy"], color=S1, lw=2.0,
            solid_capstyle="round", solid_joinstyle="round",
            label="observed occupancy", zorder=3)

    # Direct-label only the called peaks - never a number on every point.
    for _, row in (s.iterrows() if len(s) else []):
        ax.plot([row["peak_position"]], [row["occupancy"]], "o", ms=8,
                color=S1, mec=SURFACE, mew=2.0, zorder=4)
        ax.annotate(f"pos {int(row['peak_position'])}\n{row['occupancy']:.2f}",
                    (row["peak_position"], row["occupancy"]),
                    textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=8, color=INK2, zorder=5)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    ax.set_xlim(p["position"].min(), p["position"].max())
    ax.set_ylim(0, ymax)
    ax.set_xlabel("position in oligo (1-based sequence index)")
    ax.set_ylabel("occupancy")
    ax.set_title(f"{seq} - cGAS occupancy vs shuffle null", fontsize=11, loc="left")
    leg = ax.legend(frameon=False, fontsize=8, loc="upper right", labelcolor=INK2)
    for t in leg.get_texts():
        t.set_color(INK2)
    tidy(ax)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig1_occupancy_{seq}.svg")
    fig.savefig(f"{outdir}/fig1_occupancy_{seq}.png", dpi=200)
    plt.close(fig)
    p.to_csv(f"{outdir}/fig1_occupancy_{seq}.csv", index=False)


def fig_contact_map(ct: pd.DataFrame, site: str, outdir: str, vmax: float) -> None:
    """Figure 2 - residue x nucleotide contact map. Fixed scale across sites."""
    sel = ct[ct["site_id"] == site]
    if len(sel) == 0:
        return
    grid = sel.pivot_table(index="residue", columns="position",
                           values="contact_frac", aggfunc="max").sort_index()
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.grid(False)
    im = ax.imshow(grid.to_numpy(), aspect="auto", cmap=BLUES,
                   vmin=0.0, vmax=vmax, origin="lower",
                   extent=[grid.columns.min() - 0.5, grid.columns.max() + 0.5,
                           grid.index.min() - 0.5, grid.index.max() + 0.5])
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("contact fraction", color=INK2, fontsize=9)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=MUTED, length=0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=9))
    ax.set_xlabel("position in oligo (1-based sequence index)")
    ax.set_ylabel("cGAS residue index")
    ax.set_title(f"{site} - per-nucleotide contact frequency", fontsize=11, loc="left")
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig2_contacts_{site}.svg")
    fig.savefig(f"{outdir}/fig2_contacts_{site}.png", dpi=200)
    plt.close(fig)
    grid.to_csv(f"{outdir}/fig2_contacts_{site}.csv")


def fig_cooccupancy(obs: pd.DataFrame, sites: pd.DataFrame, seq: str,
                    outdir: str) -> bool:
    """Figure 5 - site x site occupancy correlation across replicates.

    Correlation has polarity (cooperative vs anti-correlated loading), so this
    is a diverging scale with a neutral midpoint, not a sequential ramp.
    """
    s = sites[sites["sequence_id"] == seq]
    if len(s) < 2:
        return False
    peaks = s["peak_position"].tolist()
    wide = (obs[(obs["sequence_id"] == seq) & (obs["position"].isin(peaks))]
            .pivot_table(index="replicate", columns="position", values="dwell_frac"))
    corr = wide.corr(method="spearman")

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    ax.grid(False)
    im = ax.imshow(corr.to_numpy(), cmap=DIVERGE,
                   norm=TwoSlopeNorm(vmin=-1, vcenter=0.0, vmax=1))
    cb = fig.colorbar(im, ax=ax, pad=0.03)
    cb.set_label("Spearman rho across replicates", color=INK2, fontsize=9)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=MUTED, length=0)
    labels = [f"pos {c}" for c in corr.columns]
    ax.set_xticks(range(len(labels)), labels, fontsize=8, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, color=INK if abs(corr.iloc[i, j]) < 0.6 else "#ffffff")
    ax.set_title(f"{seq} - site co-occupancy", fontsize=11, loc="left")
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig5_cooccupancy_{seq}.svg")
    fig.savefig(f"{outdir}/fig5_cooccupancy_{seq}.png", dpi=200)
    plt.close(fig)
    corr.to_csv(f"{outdir}/fig5_cooccupancy_{seq}.csv")
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--occupancy", required=True)
    ap.add_argument("--contacts", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    prof = pd.read_csv(args.profile, sep="\t")
    sites = pd.read_csv(args.sites, sep="\t")
    obs = pd.read_csv(args.occupancy, sep="\t")
    ct = pd.read_csv(args.contacts, sep="\t")

    if len(prof) == 0:
        raise SystemExit("ERROR: occupancy profile is empty")

    # Fixed scales across every panel.
    ymax = float(max(prof["occupancy"].max() * 1.25, prof["null_threshold"].max() * 3))
    vmax = float(ct["contact_frac"].max()) if len(ct) else 1.0

    made = []
    for seq in sorted(prof["sequence_id"].unique()):
        fig_occupancy_track(prof, sites, seq, args.outdir, ymax)
        made.append(f"fig1_occupancy_{seq}")
        if fig_cooccupancy(obs, sites, seq, args.outdir):
            made.append(f"fig5_cooccupancy_{seq}")

    for site in sorted(sites["site_id"].unique()) if len(sites) else []:
        fig_contact_map(ct, site, args.outdir, vmax)
        made.append(f"fig2_contacts_{site}")

    print(f"[make_figures] wrote {len(made)} figure set(s) to {args.outdir}")
    for m in made:
        print(f"  {m}  (.svg .png .csv)")
    if len(sites) == 0:
        print("[make_figures] no called sites: occupancy tracks written, "
              "per-site panels skipped (valid null result)")


if __name__ == "__main__":
    main()
