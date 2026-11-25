#!/usr/bin/env python3
"""Pool bound frames per called site into the pose ensembles CLUSTER_POSES reads.

Runs after site calling, for ONE (arm, sequence_id): reads every replicate's
trajectory, selects the frames in which cGAS occupies each called site, and
writes one .npz per site containing `coords`, `bead_labels` and `replicate`.

Why the DNA is in the npz
-------------------------
cluster_poses.py superposes with kabsch_superpose(coords, coords[0]) over
EVERY bead in the file and does not subset by the contacts table. If the npz
held only cGAS beads, each frame's cGAS would be aligned onto frame 0's cGAS,
which removes the very thing being clustered - where cGAS sits on the DNA - and
leaves only internal conformational change. Including the site's DNA beads
anchors the alignment, because the duplex is near-rigid and near-identical
across frames, so the residual RMSD is dominated by cGAS placement.

Consequence to record: the RMSD is computed over DNA + protein beads together,
and the DNA contributes almost nothing to it. params.rmsd_cutoff is therefore
on that combined scale, which is SMALLER than a protein-only RMSD for the same
displacement. Calibrate it from `cluster_poses.py --report-only` on these files
rather than from a protein-only intuition.

Two choices here are scientific, exposed as flags rather than buried:
  --pose-beads    which beads define a pose (protein, or interface only)
  --site-overlap  how much of a frame's footprint must fall in the site
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd

SITE_COLS = ["sequence_id", "site_id", "run_start", "run_end"]


def die(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(1)


def load_frames(path: str, n_beads: int) -> np.ndarray:
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=False)
        if "coords" not in d:
            die(f"{path} has no 'coords' array")
        return np.asarray(d["coords"], dtype=np.float64)
    try:
        import MDAnalysis as mda
    except ImportError:
        die("MDAnalysis is required to read DCD trajectories")
    u = mda.Universe.empty(n_beads, trajectory=True)
    u.load_new(path, format="DCD")
    return np.asarray([ts.positions.copy() for ts in u.trajectory],
                      dtype=np.float64)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="TSV: replicate, path. Replicate identity comes from "
                         "this column, never from the order of the files.")
    ap.add_argument("--bead-map", required=True)
    ap.add_argument("--sites", required=True, help="sites.tsv for THIS arm")
    ap.add_argument("--sequence-id", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--contact-cutoff", type=float, required=True,
                    help="Angstrom; must match the value used for occupancy, or "
                         "a frame could be 'bound' for site calling and not here")
    ap.add_argument("--footprint-min", type=int, required=True)
    ap.add_argument("--equil-frac", type=float, default=0.2)
    ap.add_argument("--site-overlap", type=float, default=0.5,
                    help="fraction of the site's positions that a frame's "
                         "footprint must cover for the frame to join that "
                         "site's ensemble (default 0.5)")
    ap.add_argument("--pose-beads", choices=["protein", "interface"],
                    default="protein",
                    help="'protein': every cGAS bead, so orientation counts. "
                         "'interface': only cGAS beads contacting the site in a "
                         "majority of its bound frames, so distal motion is "
                         "ignored. Default protein - it needs no extra "
                         "threshold and is reproducible.")
    ap.add_argument("--min-frames", type=int, default=50)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-manifest", required=True,
                    help="TSV written for the pipeline: sequence_id, site_id, npz_path")
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    if not 0.0 < args.site_overlap <= 1.0:
        die(f"--site-overlap must be in (0,1], got {args.site_overlap}")

    bm = pd.read_csv(args.bead_map, sep="\t")
    for c in ("bead_index", "role"):
        if c not in bm.columns:
            die(f"{args.bead_map} is missing column {c}")
    bm = bm.sort_values("bead_index").reset_index(drop=True)
    if list(bm["bead_index"]) != list(range(len(bm))):
        die(f"{args.bead_map} bead_index must be 0..N-1 with no gaps")

    dna = bm[bm["role"] == "dna"]
    pro = bm[bm["role"] == "protein"]
    if dna.empty or pro.empty:
        die("bead map needs both dna and protein beads")
    if dna["nt_index"].isna().any():
        die("every dna bead needs an nt_index")
    dna_nt = dna["nt_index"].astype(int).to_numpy()
    dna_idx = dna["bead_index"].to_numpy()
    pro_idx = pro["bead_index"].to_numpy()
    n_beads = len(bm)

    sites = pd.read_csv(args.sites, sep="\t")
    missing = [c for c in SITE_COLS if c not in sites.columns]
    if missing:
        die(f"{args.sites} is missing column(s): {', '.join(missing)}")
    sites = sites[sites["sequence_id"] == args.sequence_id]
    if sites.empty:
        # No called site for this sequence is a real outcome, not a failure.
        # Write the (empty) manifest so the pipeline can proceed.
        with open(args.out_manifest, "w") as fh:
            fh.write("sequence_id\tsite_id\tnpz_path\n")
        with open(args.out_summary, "w") as fh:
            json.dump({"sequence_id": args.sequence_id, "arm": args.arm,
                       "n_sites": 0, "reason": "no called site for this sequence"},
                      fh, indent=2)
        print(f"[extract_poses] {args.sequence_id} arm={args.arm}: no called sites")
        return

    man = pd.read_csv(args.manifest, sep="\t")
    for c in ("replicate", "path"):
        if c not in man.columns:
            die(f"{args.manifest} is missing column {c}")
    if man["replicate"].duplicated().any():
        die(f"{args.manifest} repeats a replicate, so frames would be labelled "
            "with the wrong replicate")

    # per-frame footprint, pooled over replicates and labelled with replicate
    frames, labels_rep = [], []
    for row in man.itertuples(index=False):
        coords = load_frames(str(row.path), n_beads)
        if coords.ndim != 3 or coords.shape[2] != 3:
            die(f"{row.path}: expected (frames, beads, 3), got {coords.shape}")
        if coords.shape[1] != n_beads:
            die(f"{row.path}: {coords.shape[1]} beads but the bead map has "
                f"{n_beads}; map and trajectory describe different systems")
        start = int(round(args.equil_frac * coords.shape[0]))
        for f in range(start, coords.shape[0]):
            frames.append(coords[f])
            labels_rep.append(int(row.replicate))
    if not frames:
        die("no post-equilibration frames across the supplied replicates")
    allf = np.asarray(frames)
    rep_arr = np.asarray(labels_rep, dtype=int)

    # which nucleotides each frame's cGAS contacts
    d = np.linalg.norm(allf[:, pro_idx, None, :] - allf[:, None, dna_idx, :],
                       axis=-1)
    close = d.min(axis=1) <= args.contact_cutoff          # (frames, dna beads)
    n_nt = int(dna_nt.max())
    nt_contact = np.zeros((allf.shape[0], n_nt + 1), dtype=bool)
    for k, nt in enumerate(dna_nt):
        nt_contact[:, nt] |= close[:, k]

    written, summaries = [], []
    for site in sites.itertuples(index=False):
        lo, hi = int(site.run_start), int(site.run_end)
        if lo > hi:
            die(f"site {site.site_id}: run_start {lo} > run_end {hi}")
        span = np.arange(max(1, lo), min(n_nt, hi) + 1)
        if span.size == 0:
            die(f"site {site.site_id}: span {lo}-{hi} lies outside 1..{n_nt}")

        covered = nt_contact[:, span].sum(axis=1)
        total = nt_contact[:, 1:].sum(axis=1)
        keep = (covered >= args.site_overlap * span.size) & \
               (total >= args.footprint_min)

        n_keep = int(keep.sum())
        if args.pose_beads == "interface":
            # Fix the bead set ONCE per site from the kept frames; it must be
            # constant across frames because bead_labels is a single array.
            if n_keep:
                hit = (d[keep][:, :, :] <= args.contact_cutoff).any(axis=2)
                frac = hit.mean(axis=0)
                sel_pro = pro_idx[frac >= 0.5]
            else:
                sel_pro = pro_idx[:0]
            if sel_pro.size == 0:
                sel_pro = pro_idx
        else:
            sel_pro = pro_idx

        site_dna = dna_idx[np.isin(dna_nt, span)]
        sel = np.concatenate([site_dna, sel_pro])
        bead_labels = np.array(
            [f"dna:{dna_nt[list(dna_idx).index(i)]}" for i in site_dna] +
            [f"pro:{i}" for i in sel_pro], dtype=object).astype("U32")

        summaries.append({"site_id": site.site_id, "n_frames_pooled": n_keep,
                          "n_replicates": int(np.unique(rep_arr[keep]).size) if n_keep else 0,
                          "n_beads": int(sel.size),
                          "below_min_frames": bool(n_keep < args.min_frames)})
        if n_keep == 0:
            print(f"[extract_poses] {site.site_id}: no frame met the site "
                  f"criterion; no npz written")
            continue

        out = f"{args.out_dir}/{site.site_id}_poses.npz"
        np.savez_compressed(out,
                            coords=allf[keep][:, sel, :].astype(np.float32),
                            bead_labels=bead_labels,
                            replicate=rep_arr[keep])
        written.append((args.sequence_id, site.site_id, out))
        print(f"[extract_poses] {site.site_id}: {n_keep} frames from "
              f"{np.unique(rep_arr[keep]).size} replicate(s), {sel.size} beads")

    with open(args.out_manifest, "w") as fh:
        fh.write("sequence_id\tsite_id\tnpz_path\n")
        for seq, sid, path in written:
            fh.write(f"{seq}\t{sid}\t{path}\n")

    with open(args.out_summary, "w") as fh:
        json.dump({"sequence_id": args.sequence_id, "arm": args.arm,
                   "pose_beads": args.pose_beads,
                   "site_overlap": args.site_overlap,
                   "frames_pooled_total": int(allf.shape[0]),
                   "n_sites": len(summaries), "sites": summaries}, fh, indent=2)


if __name__ == "__main__":
    main()
