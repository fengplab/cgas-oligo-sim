#!/usr/bin/env python3
"""Merge GENESIS-CG .gro files into one system and emit the bead map.

Two jobs that must agree with each other, so they are done from one parse:

  1. A single valid .gro. Concatenating .gro files does NOT work - each carries
     a title line, an atom count and a trailing box line, so `cat a.gro b.gro`
     produces a file with those in the middle. Atom indices must also be
     renumbered continuously, because the trajectory's bead order is this file's
     order.
  2. The bead map that extract_occupancy.py requires: bead_index (0..N-1, no
     gaps), role, nt_index, strand, copy.

Conventions below were verified against genesis_cg_tool's own reference outputs
(test/coarse-graining), not assumed:

  * Protein (AICG2+): one bead per residue, named CA. VP_cg_ref.gro has 244
    beads for 244 CA atoms.
  * DNA (3SPN.2C): beads DP/DS/DB, three per nucleotide EXCEPT the 5'-terminal
    nucleotide, which has no DP. So a strand of n nt gives 3n-1 beads and a
    duplex gives 2(3n-1). The 30-mer reference is 178 beads = 2(3*30-1).
  * Strand 2 of a duplex is written as the REVERSE COMPLEMENT in file order:
    for the 30-mer, residues 1-30 reproduce the input FASTA and residues 31-60
    are its reverse complement, both 5'->3'.

That last point drives nt_index. extract_occupancy takes n_nt = max(nt_index)
and 1-based numbering, so nt_index must be the BASE-PAIR position shared by
both strands. Strand 1 residue i is base pair i; strand 2 residue j is base
pair (n + 1 - j). Numbering strand 2 in file order instead would assign every
strand-2 bead to the mirror-image base pair and flip half the occupancy
profile end to end - silently, since the file would still be well formed.

Units: a .gro is in NANOMETRES and GENESIS multiplies by 10 on read
(molecules.fpp: `atom_coord * 10.0_wp`), while the [BOUNDARY] box_size_x/y/z in
the control file are ANGSTROM. So --box-angstrom is divided by 10 for the .gro
box line. Writing an Angstrom value straight into a .gro would inflate the
volume a thousandfold.
"""

import argparse
import sys

DNA_BEADS = {"DP", "DS", "DB"}
PROTEIN_BEADS = {"CA"}


def die(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(1)


def read_gro(path):
    """Parse a .gro the way GENESIS parses it: (i5,a5,a5,i5) then free-format.

    Mirroring fileio_grocrd.fpp matters - if this parse and GENESIS's disagree
    about a column, the bead map describes a different system from the one that
    is simulated.
    """
    lines = open(path).read().splitlines()
    if len(lines) < 3:
        die(f"{path}: too short to be a .gro")
    try:
        n = int(lines[1])
    except ValueError:
        die(f"{path}: line 2 is not an atom count: {lines[1]!r}")
    body = lines[2:2 + n]
    if len(body) != n:
        die(f"{path}: declares {n} atoms but has {len(body)} atom lines")
    box = lines[2 + n].split() if len(lines) > 2 + n else []

    atoms = []
    for ln, line in enumerate(body, start=3):
        try:
            res_no = int(line[0:5])
        except ValueError:
            die(f"{path}:{ln}: columns 1-5 are not a residue number: {line[0:5]!r}")
        res_name = line[5:10].strip()
        bead_name = line[10:15].strip()
        rest = line[20:].split()
        if len(rest) not in (3, 6):
            die(f"{path}:{ln}: expected 3 or 6 coordinate fields, got {len(rest)}")
        atoms.append({"res_no": res_no, "res_name": res_name,
                      "bead": bead_name, "xyz": rest[:3]})
    return atoms, box


def split_dna_strands(atoms, path):
    """Group DNA residues into strands.

    A strand's 5'-terminal nucleotide has no DP bead, which is a structural
    signal rather than a guess about ordering, so it is used to find strand
    starts. Residue numbering is continuous across strands in the tool's
    output, so it cannot be used for this on its own.
    """
    beads_by_res, order = {}, []
    for a in atoms:
        if a["res_no"] not in beads_by_res:
            beads_by_res[a["res_no"]] = set()
            order.append(a["res_no"])
        beads_by_res[a["res_no"]].add(a["bead"])

    strands, current = [], []
    for res_no in order:
        if "DP" not in beads_by_res[res_no]:
            if current:
                strands.append(current)
            current = [res_no]
        else:
            if not current:
                die(f"{path}: first DNA residue {res_no} has a DP bead, so no "
                    "5' terminus could be identified. The bead naming may differ "
                    "from the reference convention; inspect the .gro.")
            current.append(res_no)
    if current:
        strands.append(current)

    if len(strands) != 2:
        die(f"{path}: found {len(strands)} DNA strand(s) by 5'-terminus "
            f"detection, expected 2 for a duplex. Strand lengths: "
            f"{[len(s) for s in strands]}")
    if len(strands[0]) != len(strands[1]):
        die(f"{path}: strand lengths differ ({len(strands[0])} vs "
            f"{len(strands[1])}); this is not a blunt duplex and the base-pair "
            "mapping would be ambiguous.")
    return strands


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dna", required=True, help="DNA .gro from build_dna.jl")
    ap.add_argument("--protein", required=True, action="append",
                    help="protein .gro from aa_2_cg.jl; repeat for extra copies")
    ap.add_argument("--box-angstrom", type=float, required=True,
                    help="cubic box edge in ANGSTROM (params.cg_box_size); "
                         "divided by 10 for the .gro box line, which is nm")
    ap.add_argument("--title", default="cgas-oligo-sim merged CG system")
    ap.add_argument("--out-gro", required=True)
    ap.add_argument("--out-bead-map", required=True)
    args = ap.parse_args()

    if args.box_angstrom <= 0:
        die(f"--box-angstrom must be > 0, got {args.box_angstrom}")

    dna_atoms, _ = read_gro(args.dna)
    bad = {a["bead"] for a in dna_atoms} - DNA_BEADS
    if bad:
        die(f"{args.dna}: unexpected DNA bead name(s) {sorted(bad)}; expected "
            f"{sorted(DNA_BEADS)}")
    strands = split_dna_strands(dna_atoms, args.dna)
    n_nt = len(strands[0])

    expected = 2 * (3 * n_nt - 1)
    if len(dna_atoms) != expected:
        die(f"{args.dna}: {len(dna_atoms)} DNA beads but a {n_nt}-mer duplex "
            f"should give 2*(3*{n_nt}-1) = {expected}. The 5'-no-DP convention "
            "may not hold for this file.")

    # base-pair position per residue; strand 2 runs antiparallel
    bp_of_res, strand_of_res = {}, {}
    for i, res_no in enumerate(strands[0], start=1):
        bp_of_res[res_no] = i
        strand_of_res[res_no] = 1
    for j, res_no in enumerate(strands[1], start=1):
        bp_of_res[res_no] = n_nt + 1 - j
        strand_of_res[res_no] = 2

    rows, merged, bead_index, res_offset = [], [], 0, 0

    for a in dna_atoms:
        merged.append((a["res_no"], a["res_name"], a["bead"], a["xyz"]))
        rows.append((bead_index, "dna", bp_of_res[a["res_no"]],
                     strand_of_res[a["res_no"]], 1))
        bead_index += 1
    res_offset = max(a["res_no"] for a in dna_atoms)

    for copy_no, path in enumerate(args.protein, start=1):
        patoms, _ = read_gro(path)
        bad = {a["bead"] for a in patoms} - PROTEIN_BEADS
        if bad:
            die(f"{path}: unexpected protein bead name(s) {sorted(bad)}; "
                f"expected {sorted(PROTEIN_BEADS)}. AICG2+ is one CA bead per "
                "residue; a side-chain model would need a different bead map.")
        n_res = len({a["res_no"] for a in patoms})
        if n_res != len(patoms):
            die(f"{path}: {len(patoms)} beads over {n_res} residues; AICG2+ "
                "must be 1:1")
        for a in patoms:
            merged.append((a["res_no"] + res_offset, a["res_name"],
                           a["bead"], a["xyz"]))
            rows.append((bead_index, "protein", "", "", copy_no))
            bead_index += 1
        res_offset += max(a["res_no"] for a in patoms)

    box_nm = args.box_angstrom / 10.0
    with open(args.out_gro, "w") as fh:
        fh.write(f"{args.title}\n")
        fh.write(f"{len(merged):>12d}\n")
        for i, (res_no, res_name, bead, xyz) in enumerate(merged, start=1):
            # i5 fields wrap at 100000 in the GRO format; GENESIS reads them
            # with (i5,a5,a5,i5), so an overflow would silently corrupt indices.
            if res_no >= 100000 or i >= 100000:
                die("residue or atom index exceeds the 5-column GRO field")
            fh.write(f"{res_no:>5d}{res_name:>5s}{bead:>5s}{i:>5d}"
                     f"{float(xyz[0]):>8.3f}{float(xyz[1]):>8.3f}"
                     f"{float(xyz[2]):>8.3f}\n")
        fh.write(f"{box_nm:>10.5f}{box_nm:>10.5f}{box_nm:>10.5f}\n")

    with open(args.out_bead_map, "w") as fh:
        fh.write("bead_index\trole\tnt_index\tstrand\tcopy\n")
        for bi, role, nt, strand, copy_no in rows:
            fh.write(f"{bi}\t{role}\t{nt}\t{strand}\t{copy_no}\n")

    n_dna = sum(1 for r in rows if r[1] == "dna")
    n_pro = sum(1 for r in rows if r[1] == "protein")
    print(f"beads in : dna {n_dna}, protein {n_pro} "
          f"({len(args.protein)} copy/copies)")
    print(f"beads out: {len(merged)}  (bead_index 0..{len(merged) - 1})")
    print(f"duplex   : {n_nt} bp; nt_index 1..{n_nt}, strand 2 reversed")
    print(f"box      : {args.box_angstrom:g} A written as {box_nm:g} nm")
    print(f"wrote {args.out_gro} and {args.out_bead_map}")


if __name__ == "__main__":
    main()
