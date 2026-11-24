#!/usr/bin/env python3
"""Synthetic protein mmCIF carrying KNOWN engineered substitutions.

Reference is 40 residues. The structure differs at author positions 187 and
195 (offset so the numbering is not 1..40), mimicking a crystallisation
construct: LYS->ASN and LEU->ARG. prepare_structures.py must find exactly
those two by alignment, with no mutation table given to it.
"""
import gemmi

REF = "MQPWHGKAMQRASEILLARKKRPVSRHIKRLCYMATQLLRD"  # 41 aa
MUT_AT = {6: "N", 15: "R"}  # 0-based into REF: K->N, L->R
AUTH_START = 180
ATOMS = {"N": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
         "R": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
         "G": ["N", "CA", "C", "O"]}
DEFAULT = ["N", "CA", "C", "O", "CB", "CG"]
ONE3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
        "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
        "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

obs = list(REF)
for i, aa in MUT_AT.items():
    obs[i] = aa

st = gemmi.Structure()
st.name = "FIXTURE"
st.spacegroup_hm = "P 1"
st.cell = gemmi.UnitCell(200, 200, 200, 90, 90, 90)
model = gemmi.Model("1")
chain = gemmi.Chain("A")
for i, aa in enumerate(obs):
    r = gemmi.Residue()
    r.name = ONE3[aa]
    r.seqid = gemmi.SeqId(AUTH_START + i, " ")
    r.het_flag = "A"
    for j, an in enumerate(ATOMS.get(aa, DEFAULT)):
        a = gemmi.Atom()
        a.name = an
        a.element = gemmi.Element("N" if an.startswith("N") else
                                  "O" if an.startswith("O") else "C")
        a.pos = gemmi.Position(i * 3.8 + j * 0.4, j * 0.5, 0.0)
        a.occ, a.b_iso = 1.0, 20.0
        r.add_atom(a)
    chain.add_residue(r)
model.add_chain(chain)
st.add_model(model)
st.setup_entities()
st.make_mmcif_document().write_file("tests/fixtures/structure_fixture.cif")

with open("tests/fixtures/structure_reference.fasta", "w") as fh:
    fh.write(">FIXTURE_REF synthetic\n")
    for i in range(0, len(REF), 60):
        fh.write(REF[i:i+60] + "\n")

with open("tests/fixtures/structures.tsv", "w") as fh:
    fh.write("structure_id\tpdb_id\tchain\tcif_path\trevert\n")
    fh.write("FIXWT\tFIXT\tA\ttests/fixtures/structure_fixture.cif\tauto\n")

exp = ";".join(f"{MUT_AT[i]}{AUTH_START+i}{REF[i]}" for i in sorted(MUT_AT))
print(f"fixture written; expected reversions: {exp}")
