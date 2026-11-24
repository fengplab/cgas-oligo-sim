#!/usr/bin/env python3
"""A structure whose chain NAMES repeat, as gemmi presents a real mmCIF.

gemmi splits one mmCIF auth chain into several Chain objects - polymer, het,
water - that all carry the same name. Real 6CT9 therefore yields two chains
called "A": the 348-residue protein and a 94-residue het/water part. After
remove_ligands_and_waters() the second is empty but still present, and
remove_chain() deletes the FIRST match by name, so dropping the empty "A"
deleted the protein.

Writing this structure to mmCIF and reading it back MERGES the parts, so the
fixture must be built in memory and handed to the code directly. That is
exactly why no earlier fixture caught the bug.
"""
import gemmi

ONE3 = {"M":"MET","Q":"GLN","P":"PRO","W":"TRP","H":"HIS","G":"GLY","K":"LYS",
        "A":"ALA","R":"ARG","S":"SER","E":"GLU","T":"THR"}
PROT = "MQPWHGKAMQRASEAGATAPK"


def _atom(name, elem, x, y=0.0):
    a = gemmi.Atom()
    a.name, a.element = name, gemmi.Element(elem)
    a.pos = gemmi.Position(x, y, 0.0)
    a.occ, a.b_iso = 1.0, 20.0
    return a


def build() -> gemmi.Structure:
    st = gemmi.Structure()
    st.name = "DUPCHAIN"
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(200, 200, 200, 90, 90, 90)
    model = gemmi.Model("1")

    ch = gemmi.Chain("A")
    for i, aa in enumerate(PROT):
        r = gemmi.Residue()
        r.name, r.seqid, r.het_flag = ONE3[aa], gemmi.SeqId(157 + i, " "), "A"
        for j, an in enumerate(["N", "CA", "C", "O", "CB"]):
            r.add_atom(_atom(an, "C", i * 3.8 + j * 0.3, j * 0.4))
        ch.add_residue(r)
    model.add_chain(ch)

    for cn, bases in (("B", "CGTCTTCGGCAAT"), ("C", "AAATTGCC")):
        ch = gemmi.Chain(cn)
        for i, b in enumerate(bases):
            r = gemmi.Residue()
            r.name, r.seqid, r.het_flag = "D" + b, gemmi.SeqId(i + 1, " "), "A"
            r.add_atom(_atom("P", "P", i * 3.4, 20.0))
            ch.add_residue(r)
        model.add_chain(ch)

    # The het and water parts, sharing names with the polymer chains.
    for cn in ("A", "B", "C"):
        ch = gemmi.Chain(cn)
        r = gemmi.Residue()
        r.name, r.seqid, r.het_flag = "HOH", gemmi.SeqId(701, " "), "H"
        r.add_atom(_atom("O", "O", 50.0))
        ch.add_residue(r)
        model.add_chain(ch)

    st.add_model(model)
    return st


if __name__ == "__main__":
    st = build()
    print("  chain objects:", [(c.name, len(c)) for c in st[0]])
    names = [c.name for c in st[0]]
    assert len(names) != len(set(names)), "fixture must have duplicate names"
    print("  duplicate names present, as intended")
