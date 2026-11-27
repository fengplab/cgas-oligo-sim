#!/usr/bin/env python3
"""Write a small synthetic mmCIF carrying deliberate unmodelled runs.

Nothing here is real cGAS: it is a 20-residue peptide plus a short duplex, with
author numbering offset from SEQRES numbering so that a script which confuses
the two produces visibly wrong answers. Used by tests/run_tests.sh.

Gaps planted (author numbering, offset = +100 from SEQRES index):
    101        1 residue   N-terminal   -> ineligible, no N-side anchor
    106-108    3 residues  internal     -> eligible
    112-113    2 residues  internal     -> eligible
    120        1 residue   C-terminal   -> ineligible, no C-side anchor
"""
import math
import sys
from pathlib import Path

SEQ = "AGKLMSTVWYACDEFGHIKL"          # 20 residues, SEQRES index 1..20
OFFSET = 100                            # auth_seq = index + OFFSET
UNOBSERVED = {1, 6, 7, 8, 12, 13, 20}   # SEQRES indices with no coordinates

THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
         "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
         "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
         "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}

DNA = "ATGC"                            # one short strand, chain B


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "fixture.cif")
    L: list[str] = []
    L.append("data_FAKE")
    L.append("#")
    L.append("_entry.id FAKE")
    L.append("#")

    # Two entities: 1 = peptide (chain A), 2 = DNA (chain B). find_gaps must
    # report gaps for the peptide only.
    L.append("loop_")
    L.append("_entity.id")
    L.append("_entity.type")
    L.append("1 polymer")
    L.append("2 polymer")
    L.append("#")
    L.append("loop_")
    L.append("_entity_poly.entity_id")
    L.append("_entity_poly.type")
    L.append("1 'polypeptide(L)'")
    L.append("2 'polydeoxyribonucleotide'")
    L.append("#")

    # The deposited SEQRES-to-coordinate mapping. pdb_mon_id is '?' exactly for
    # residues with no coordinates, which is what find_gaps reads.
    L.append("loop_")
    for tag in ("asym_id", "entity_id", "seq_id", "mon_id", "ndb_seq_num",
                "pdb_seq_num", "auth_seq_num", "pdb_mon_id", "auth_mon_id",
                "pdb_strand_id", "pdb_ins_code", "hetero"):
        L.append(f"_pdbx_poly_seq_scheme.{tag}")
    for i, aa in enumerate(SEQ, start=1):
        three = THREE[aa]
        auth = i + OFFSET
        seen = i not in UNOBSERVED
        mon = three if seen else "?"
        L.append(f"A 1 {i} {three} {i} {auth} {auth} {mon} {mon} A . n")
    for j, nt in enumerate(DNA, start=1):
        L.append(f"B 2 {j} D{nt} {j} {j} {j} D{nt} D{nt} B . n")
    L.append("#")

    # Two categories gemmi does not model and a PDB round trip would destroy.
    # promote.py must preserve them: prepare_structures.py cross-checks its WT
    # reversion against _struct_ref_seq_dif, and the dimer row expands
    # _pdbx_struct_assembly. The test asserts both survive.
    L.append("loop_")
    for tag in ("align_id", "pdbx_pdb_id_code", "mon_id", "pdbx_pdb_strand_id",
                "seq_num", "pdbx_auth_seq_num", "db_mon_id", "pdbx_seq_db_name",
                "pdbx_seq_db_accession_code", "details"):
        L.append(f"_struct_ref_seq_dif.{tag}")
    # SEQRES position 4 is LEU here while the reference database says MET: a
    # stand-in for an engineered mutation.
    # prepare_structures.py's cross-check reads strand, auth_seq, db_mon_id,
    # mon_id AND details together, and gemmi returns no rows at all if one tag
    # is missing - the check then reports None instead of failing, so the
    # fixture carries every column a deposited file would.
    L.append("1 FAKE LEU A 4 104 MET UNP FAKE01 'engineered mutation'")
    L.append("#")
    L.append("loop_")
    L.append("_pdbx_struct_assembly.id")
    L.append("_pdbx_struct_assembly.details")
    L.append("_pdbx_struct_assembly.oligomeric_details")
    L.append("_pdbx_struct_assembly.oligomeric_count")
    L.append("1 author_defined_assembly dimeric 2")
    L.append("#")

    # Backbone-only coordinates for the observed peptide residues, and one atom
    # per nucleotide. Geometry is a straight line: no script under test reads
    # it, and a bogus bond length is less misleading than a plausible one.
    L.append("loop_")
    # Tag order and membership follow what gemmi itself writes: it reads a
    # fixed set, and omitting label_alt_id or pdbx_PDB_ins_code yields a
    # document that parses as a table but produces ZERO models.
    for tag in ("group_PDB", "id", "type_symbol", "label_atom_id",
                "label_alt_id", "label_comp_id", "label_asym_id",
                "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
                "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
                "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id",
                "auth_asym_id", "pdbx_PDB_model_num"):
        L.append(f"_atom_site.{tag}")

    # One CA every 3.80 A along a slowly turning path. Exactly 3.80 matters:
    # make_fake_models.py interpolates inserted residues across a gap and the
    # resulting C-N distances follow directly from CA spacing. A straight line
    # would also give 3.80, but a collinear point set makes Kabsch
    # superposition degenerate, so the path turns.
    ca = []
    pos = [0.0, 0.0, 0.0]
    for i in range(1, len(SEQ) + 1):
        t = 0.25 * i
        u = (math.cos(t) * 0.6 + 0.8, math.sin(t) * 0.6, math.sin(t * 0.7) * 0.3)
        n = math.sqrt(sum(c * c for c in u))
        pos = [pos[k] + 3.80 * u[k] / n for k in range(3)]
        ca.append(tuple(pos))

    def backbone(i):
        """N, CA, C, O for residue i (1-based), oriented along the local path."""
        p = ca[i - 1]
        nxt = ca[i] if i < len(ca) else ca[i - 2]
        d = [nxt[k] - p[k] for k in range(3)]
        n = math.sqrt(sum(c * c for c in d)) or 1.0
        u = [c / n for c in d]
        if i == len(ca):
            u = [-c for c in u]
        return {
            "N": [p[k] - 0.95 * u[k] for k in range(3)],
            "CA": list(p),
            "C": [p[k] + 1.52 * u[k] for k in range(3)],
            "O": [p[k] + 2.10 * u[k] for k in range(3)],
        }

    serial = 0
    for i, aa in enumerate(SEQ, start=1):
        if i in UNOBSERVED:
            continue
        three, auth = THREE[aa], i + OFFSET
        bb = backbone(i)
        for name, elem in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
            serial += 1
            x, y, z = bb[name]
            L.append(f"ATOM {serial} {elem} {name} . {three} A 1 {i} ? "
                     f"{x:.3f} {y:.3f} {z:.3f} 1.00 20.00 ? {auth} A 1")
    for j, nt in enumerate(DNA, start=1):
        serial += 1
        L.append(f"ATOM {serial} P P . D{nt} B 2 {j} ? "
                 f"{j * 6.0:.3f} 20.000 0.000 1.00 25.00 ? {j} B 1")
    L.append("#")

    out.write_text("\n".join(L) + "\n")

    # Ground truth for every residue, observed or not. make_fake_models.py uses
    # it so the inserted residues sit on the same 3.80 A path as the rest and
    # the peptide bonds come out right by construction; without it the fixture
    # has to interpolate across a chord that is shorter than the path it
    # replaces, which compresses the bonds it is meant to be testing.
    truth = out.with_suffix(".truth.tsv")
    rows = ["auth_seq\tx\ty\tz"]
    for i in range(1, len(SEQ) + 1):
        x, y, z = ca[i - 1]
        rows.append(f"{i + OFFSET}\t{x:.3f}\t{y:.3f}\t{z:.3f}")
    truth.write_text("\n".join(rows) + "\n")
    print(f"wrote {out} and {truth}")


if __name__ == "__main__":
    main()
