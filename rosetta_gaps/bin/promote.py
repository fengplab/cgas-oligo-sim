#!/usr/bin/env python3
"""Write a chosen Rosetta model back into the deposited mmCIF, for structures.tsv.

The repaired file is the ORIGINAL mmCIF with the protein chain's coordinates
replaced, not a fresh file built from the model. That matters because
bin/prepare_structures.py reads categories a PDB-to-mmCIF conversion would
throw away:

  _struct_ref_seq_dif     the deposited substitution records, which the WT
                          reversion cross-checks its alignment against. Lose
                          them and crosscheck_agrees silently becomes None.
  _pdbx_struct_assembly*  what the dimer row expands with transform_to_assembly.
                          Lose them and that row dies.
  _pdbx_poly_seq_scheme   the SEQRES-to-coordinate map. Left stale it would
                          still describe the repaired residues as unobserved,
                          and find_gaps.py would offer to model them again.
                          This script updates it.

Everything except the protein chain is taken from the source file untouched, so
the DNA and the zinc are the deposited coordinates and not Rosetta's view of
them. Side chains NEAR a rebuilt loop do come from Rosetta: loopmodel repacks
loop neighbours unless -loops:fix_natsc is set, and it is not set here. For a
C-alpha coarse-grained model that is immaterial, but it is not the crystal.

Promotion is deliberately a separate command. One model becomes the native
state of the AICG2+ reference and nothing downstream can penalise it, so the
choice stays an explicit act with the spread report in front of you.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import gemmi
import pandas as pd


def die(msg: str) -> None:
    sys.exit(f"promote: {msg}")


def read_model(path: Path) -> gemmi.Structure:
    if path.suffix == ".gz":
        st = gemmi.read_pdb_string(gzip.decompress(path.read_bytes()).decode())
    else:
        st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--model", required=True,
                    help="the model to promote, from models/ (.pdb or .pdb.gz)")
    ap.add_argument("--out", required=True,
                    help="output mmCIF, e.g. data/structures/6ct9.repaired.cif")
    ap.add_argument("--force", action="store_true",
                    help="write even though the spread report calls a rebuilt "
                         "gap undetermined")
    args = ap.parse_args()

    work = Path(args.work_dir)
    prov = json.loads((work / "provenance.json").read_text())
    numbering = pd.read_csv(work / "numbering.tsv", sep="\t",
                            dtype=str).fillna("")
    expected = (work / "expected_sequence.txt").read_text().strip()
    chain_name = prov["chain"]

    report = work / "ngk_report.tsv"
    if report.exists():
        rep = pd.read_csv(report, sep="\t", dtype=str).fillna("")
        # check_models records each model by its path relative to the work
        # directory, because two array tasks both produce a "..._0001". Accept
        # either that label or a bare filename, but only when the filename is
        # unambiguous - promoting the wrong task's model would be invisible.
        target = Path(args.model)
        try:
            rel = str(target.resolve().relative_to(work.resolve()))
        except ValueError:
            rel = target.name
        row = rep[rep["model"] == rel]
        if row.empty:
            by_name = rep[rep["model"].apply(lambda m: Path(m).name == target.name)]
            if len(by_name) > 1:
                die(f"{target.name} matches {len(by_name)} rows in {report} "
                    f"({', '.join(by_name['model'])}). Pass the path as it "
                    "appears in the report.")
            row = by_name
        if row.empty:
            die(f"{Path(args.model).name} is not in {report}. Run "
                "check_models.py --stage ngk first: an unchecked model may "
                "carry a chain break or a shifted fixed region.")
        if row.iloc[0]["valid"] != "yes":
            die(f"{Path(args.model).name} FAILED validation: "
                f"{row.iloc[0]['problems']}")
    else:
        die(f"no {report}. Run check_models.py --stage ngk before promoting.")

    spread = work / "spread_report.tsv"
    if spread.exists():
        sp = pd.read_csv(spread, sep="\t", dtype=str).fillna("")
        undetermined = sp[sp["verdict"].str.contains("NOT determined")]
        if len(undetermined) and not args.force:
            print(undetermined.to_string(index=False), file=sys.stderr)
            die("the gap(s) above are not determined by the data: the top "
                "models disagree, so promoting one asserts a conformation the "
                "structure does not support. Decide deliberately and pass "
                "--force, and record it as a limitation.")

    model = read_model(Path(args.model))
    src_doc = gemmi.cif.read(prov["source_cif"])
    block = src_doc.sole_block()
    src = gemmi.make_structure_from_block(block)
    src.setup_entities()

    mchain = None
    for c in model[0]:
        if c.name == chain_name:
            mchain = c
    if mchain is None:
        die(f"chain {chain_name} absent from {args.model}")
    if len(mchain) != len(numbering):
        die(f"{args.model} chain {chain_name} has {len(mchain)} residues, "
            f"numbering.tsv has {len(numbering)}. Refusing to map them.")

    # Rebuild the protein chain in author numbering, from the model.
    new_chain = gemmi.Chain(chain_name)
    rows = list(numbering.itertuples(index=False))
    for res, nrow in zip(mchain, rows):
        copy = gemmi.Residue()
        copy.name = res.name
        copy.seqid = gemmi.SeqId(str(nrow.auth_seq) +
                                 (nrow.ins_code if nrow.ins_code else ""))
        copy.het_flag = "A"
        for atom in res:
            copy.add_atom(atom)
        new_chain.add_residue(copy)

    seq_in = "".join(
        (lambda i: i.one_letter_code.upper() if i and i.one_letter_code else "X")
        (gemmi.find_tabulated_residue(r.name)) for r in new_chain)
    if seq_in != expected:
        die(f"the promoted chain's sequence does not match "
            f"expected_sequence.txt ({len(seq_in)} vs {len(expected)}).")

    model_idx = [i for i, c in enumerate(src[0]) if c.name == chain_name]
    if len(model_idx) != 1:
        die(f"source file has {len(model_idx)} chains named {chain_name}; "
            "expected exactly one.")
    src[0].add_chain(new_chain, pos=model_idx[0])
    del src[0][model_idx[0] + 1]
    src.setup_entities()

    # Coordinates only: leaving `assembly` off keeps the deposited assembly,
    # substitution and SEQRES categories exactly as they were.
    groups = gemmi.MmcifOutputGroups(False)
    groups.atoms = True
    groups.atom_type = True
    src.update_mmcif_block(block, groups)

    # The scheme must stop calling the repaired residues unobserved.
    table = block.find("_pdbx_poly_seq_scheme.",
                       ["pdb_strand_id", "pdb_seq_num", "mon_id",
                        "pdb_mon_id", "auth_mon_id"])
    if len(table):
        want = {str(r.auth_seq): r for r in rows if r.modelled == "1"}
        touched = 0
        for r in table:
            if r.str(0) != chain_name:
                continue
            if r.str(1) in want and r.str(3) in ("?", ".", ""):
                r[3] = r.str(2)
                r[4] = r.str(2)
                touched += 1
        print(f"_pdbx_poly_seq_scheme: marked {touched} residue(s) observed")
    else:
        print("warning: source file has no _pdbx_poly_seq_scheme to update; "
              "find_gaps.py will refuse to read the repaired file",
              file=sys.stderr)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    src_doc.write_file(str(out))

    # Read it back: a file that cannot be re-read is not a deliverable.
    check = gemmi.cif.read(str(out))
    cblock = check.sole_block()
    cst = gemmi.make_structure_from_block(cblock)
    cst.setup_entities()
    got = None
    for c in cst[0]:
        if c.name == chain_name:
            got = c
    if got is None or len(got) != len(numbering):
        die(f"{out} re-read with "
            f"{0 if got is None else len(got)} residues in chain "
            f"{chain_name}, expected {len(numbering)}.")
    for cat, why in (("_struct_ref_seq_dif.", "WT-reversion cross-check"),
                     ("_pdbx_struct_assembly.", "assembly expansion")):
        now = len(cblock.find_mmcif_category(cat))
        if now == 0:
            print(f"note: {out.name} carries no {cat} - the source file had "
                  f"none either, so {why} has nothing to use.")

    print(f"\npromoted {Path(args.model).name} -> {out}")
    print("structures.tsv line (check revert and assembly against the row you "
          "are replacing):")
    print(f"\n{prov['structure_id']}\t{prov.get('pdb_id','')}\t{chain_name}\t"
          f"{out}\tauto\t")
    print("\nThe WT reversion has NOT been applied: run prepare_structures.py "
          "on this file as usual. That order is deliberate - reverting before "
          "the repair would leave any substitution the rebuild introduced "
          "undetected.")


if __name__ == "__main__":
    main()
