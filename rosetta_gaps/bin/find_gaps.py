#!/usr/bin/env python3
"""Enumerate unmodelled residue runs in the structures named in structures.tsv.

Reads the deposited SEQRES-to-coordinate mapping (_pdbx_poly_seq_scheme) rather
than inferring numbering, and cross-checks it against the residues that
actually carry atoms. A gap's residue identities and author numbering both come
from the deposition, so nothing here is transcribed by hand.

Writes one row per gap: author numbering (what you read in the PDB), the
missing one-letter sequence, the flanking observed residues, whether the gap is
internal or terminal, and whether it is eligible for KIC loop modelling.

Coordinate conventions: mmCIF author numbering is 1-based and may be
non-contiguous; _pdbx_poly_seq_scheme.seq_id is the 1-based SEQRES index of the
construct. Neither is Rosetta pose numbering - make_inputs.py does that
conversion and emits its own map.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gemmi
import pandas as pd

# Terminal runs cannot be built by KIC: the C-alpha atoms of the first and last
# residue of a loop definition stay fixed, so a terminal run has no anchor on
# one side. Long runs are excluded because closure quality falls off with
# length; the cutoff is a parameter, not a law.
DEFAULT_MAX_LEN = 12


def die(msg: str) -> None:
    sys.exit(f"find_gaps: {msg}")


def read_sheet(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    df.columns = [c.strip() for c in df.columns]
    for col in ("structure_id", "chain", "cif_path"):
        if col not in df.columns:
            die(f"{path}: missing required column {col!r}")
    if df.empty:
        die(f"{path}: no uncommented rows. Nothing to do.")
    return df


def poly_seq_scheme(doc: gemmi.cif.Document, cif: str) -> list[dict]:
    """Every SEQRES position of every polymer chain, observed or not."""
    block = doc.sole_block()
    cols = ["pdb_strand_id", "entity_id", "seq_id", "mon_id",
            "pdb_seq_num", "pdb_mon_id", "pdb_ins_code"]
    table = block.find("_pdbx_poly_seq_scheme.", cols)
    if not len(table):
        die(f"{cif}: no _pdbx_poly_seq_scheme category. This script will not "
            "guess author numbering for unobserved residues from the "
            "flanking offsets - an insertion code or a renumbered deposit "
            "would make that guess wrong without any error. Re-download the "
            "mmCIF from https://files.rcsb.org/download/<ID>.cif, which "
            "carries the category.")
    rows = []
    for r in table:
        rows.append({
            "chain": r.str(0),
            "entity": r.str(1),
            "seq_id": int(r.str(2)),
            # mon_id is the SEQRES identity and is always present; pdb_mon_id
            # is '?' exactly when the residue has no coordinates.
            "mon_id": r.str(3),
            "auth_seq": r.str(4),
            "observed": r.str(5) not in ("?", ".", ""),
            "ins_code": "" if r.str(6) in ("?", ".") else r.str(6),
        })
    return rows


def residues_with_atoms(st: gemmi.Structure) -> set[tuple[str, str]]:
    """(chain, auth_seq) for every residue that actually carries atoms."""
    present = set()
    for chain in st[0]:
        for res in chain:
            if len(res):
                present.add((chain.name, str(res.seqid.num)))
    return present


def protein_chains(st: gemmi.Structure) -> set[str]:
    """Author chain names whose polymer is peptide."""
    names = set()
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) and poly.check_polymer_type() in (
            gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD,
        ):
            names.add(chain.name)
    return names


def one_letter(mon_id: str) -> str:
    code = gemmi.find_tabulated_residue(mon_id)
    letter = code.one_letter_code if code else "?"
    # gemmi lowercases the code for non-standard residues; a gap containing a
    # modified residue is worth seeing rather than silently uppercasing.
    return letter.upper() if letter and letter.isalpha() else "X"


def runs(scheme: list[dict], chain: str) -> list[dict]:
    """Contiguous unobserved runs within one chain, in SEQRES order."""
    rows = sorted([r for r in scheme if r["chain"] == chain],
                  key=lambda r: r["seq_id"])
    if not rows:
        return []
    out, current = [], None
    for r in rows:
        if not r["observed"]:
            if current is None:
                current = {"start_idx": r["seq_id"], "end_idx": r["seq_id"],
                           "mon": [r["mon_id"]], "auth": [r["auth_seq"]],
                           "ins": [r["ins_code"]]}
            else:
                current["end_idx"] = r["seq_id"]
                current["mon"].append(r["mon_id"])
                current["auth"].append(r["auth_seq"])
                current["ins"].append(r["ins_code"])
        elif current is not None:
            out.append(current)
            current = None
    if current is not None:
        out.append(current)

    first_obs = min((r["seq_id"] for r in rows if r["observed"]), default=None)
    last_obs = max((r["seq_id"] for r in rows if r["observed"]), default=None)
    if first_obs is None:
        die(f"chain {chain}: no observed residues at all.")

    by_idx = {r["seq_id"]: r for r in rows}
    for run in out:
        run["position"] = ("n_term" if run["end_idx"] < first_obs else
                           "c_term" if run["start_idx"] > last_obs else
                           "internal")
        n = by_idx.get(run["start_idx"] - 1)
        c = by_idx.get(run["end_idx"] + 1)
        run["flank_n"] = n["auth_seq"] if n and n["observed"] else ""
        run["flank_c"] = c["auth_seq"] if c and c["observed"] else ""
        run["flank_n_aa"] = one_letter(n["mon_id"]) if n else ""
        run["flank_c_aa"] = one_letter(c["mon_id"]) if c else ""
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True,
                    help="structures.tsv. Commented rows are skipped, so this "
                         "reads exactly the structures you have uncommented.")
    ap.add_argument("--out", required=True, help="gaps TSV to write")
    ap.add_argument("--max-len", type=int, default=DEFAULT_MAX_LEN,
                    help=f"longest run to mark eligible (default {DEFAULT_MAX_LEN})")
    args = ap.parse_args()

    sheet = read_sheet(Path(args.samplesheet))
    rows: list[dict] = []

    for _, row in sheet.iterrows():
        sid, cif = row["structure_id"], row["cif_path"]
        if not cif:
            die(f"{sid}: cif_path is empty.")
        if not Path(cif).exists():
            die(f"{sid}: cif_path does not exist: {cif}")

        doc = gemmi.cif.read(cif)
        st = gemmi.make_structure_from_block(doc.sole_block())
        st.setup_entities()
        if not len(st):
            die(f"{sid}: {cif} parsed to zero models.")

        scheme = poly_seq_scheme(doc, cif)
        prot = protein_chains(st)
        if not prot:
            die(f"{sid}: no peptide chain in {cif}.")

        wanted = (sorted(prot) if row["chain"].strip().lower() == "all"
                  else [c.strip() for c in row["chain"].split(",") if c.strip()])
        unknown = [c for c in wanted if c not in prot]
        if unknown:
            die(f"{sid}: chain(s) {unknown} are not peptide chains in {cif}. "
                f"Peptide chains present: {sorted(prot)}")

        with_atoms = residues_with_atoms(st)

        for chain in wanted:
            # The scheme says which residues were deposited with coordinates;
            # confirm against the model. A residue the scheme calls observed
            # but that carries no atoms would otherwise be modelled as present
            # and then fail inside Rosetta with no explanation.
            claimed = {(r["chain"], r["auth_seq"]) for r in scheme
                       if r["chain"] == chain and r["observed"]}
            missing_atoms = claimed - with_atoms
            if missing_atoms:
                die(f"{sid} chain {chain}: _pdbx_poly_seq_scheme lists "
                    f"{len(missing_atoms)} residue(s) as observed that carry "
                    f"no atoms in the model, e.g. "
                    f"{sorted(missing_atoms)[:5]}. The file and its own "
                    "mapping disagree; do not model against it.")

            for run in runs(scheme, chain):
                length = run["end_idx"] - run["start_idx"] + 1
                if any(i for i in run["ins"] if i):
                    eligible, reason = "no", "insertion code inside the run"
                elif run["position"] != "internal":
                    eligible, reason = "no", (
                        f"{run['position']}: KIC holds the C-alpha of the "
                        "first and last loop residue fixed, so a terminal run "
                        "has no anchor on one side")
                elif length > args.max_len:
                    eligible, reason = "no", (
                        f"{length} residues exceeds --max-len {args.max_len}")
                elif not run["flank_n"] or not run["flank_c"]:
                    eligible, reason = "no", "a flanking residue is unobserved"
                else:
                    eligible, reason = "yes", ""

                rows.append({
                    "structure_id": sid,
                    "pdb_id": row.get("pdb_id", ""),
                    "chain": chain,
                    "gap_start_auth": run["auth"][0],
                    "gap_end_auth": run["auth"][-1],
                    "length": length,
                    "missing_seq": "".join(one_letter(m) for m in run["mon"]),
                    "seqres_start": run["start_idx"],
                    "seqres_end": run["end_idx"],
                    "flank_n_auth": run["flank_n"],
                    "flank_n_aa": run["flank_n_aa"],
                    "flank_c_auth": run["flank_c"],
                    "flank_c_aa": run["flank_c_aa"],
                    "position": run["position"],
                    "eligible": eligible,
                    "reason": reason,
                    "cif_path": cif,
                })

    if not rows:
        print("find_gaps: no unmodelled runs found in any listed chain. "
              "Nothing to model.", file=sys.stderr)

    cols = ["structure_id", "pdb_id", "chain", "gap_start_auth", "gap_end_auth",
            "length", "missing_seq", "seqres_start", "seqres_end",
            "flank_n_auth", "flank_n_aa", "flank_c_auth", "flank_c_aa",
            "position", "eligible", "reason", "cif_path"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=cols).to_csv(out, sep="\t", index=False)

    n_elig = sum(1 for r in rows if r["eligible"] == "yes")
    print(f"find_gaps: {len(rows)} run(s) across "
          f"{len({(r['structure_id'], r['chain']) for r in rows})} chain(s); "
          f"{n_elig} eligible. -> {out}")
    for r in rows:
        mark = "  " if r["eligible"] == "yes" else "x "
        print(f"{mark}{r['structure_id']:<22} {r['chain']:<3} "
              f"{r['gap_start_auth']}-{r['gap_end_auth']} "
              f"len {r['length']:<3} {r['missing_seq']:<20} {r['reason']}")


if __name__ == "__main__":
    main()
