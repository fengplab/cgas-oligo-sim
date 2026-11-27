#!/bin/bash
# Everything in rosetta_gaps/ that can be tested without Rosetta.
#
#     bash rosetta_gaps/tests/run_tests.sh [scratch dir]
#
# Covers: gap enumeration and its four refusal paths, job generation and its
# numbering, the three model-validity checks, the spread arithmetic, promotion
# with metadata preservation, and acceptance of the promoted file by
# bin/prepare_structures.py.
#
# NOT covered, because Rosetta is not here: remodel, loopmodel, the flag names
# this repository writes into flag files, container paths, and Slurm.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
RG=$(dirname "$HERE")
REPO=$(dirname "$RG")
T="${1:-$(mktemp -d)}"
mkdir -p "$T"
PASS=0
FAIL=0

ok()   { PASS=$((PASS+1)); echo "ok   $1"; }
bad()  { FAIL=$((FAIL+1)); echo "FAIL $1"; }
is()   { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (got '$2', want '$3')"; fi; }
has()  { if grep -q -- "$2" "$3" 2>/dev/null; then ok "$1"; else bad "$1"; fi; }
dies() { if "$@" >/dev/null 2>&1; then bad "should have failed: $*"; else ok "refused: ${*:1:3}"; fi; }

echo "scratch: $T"
echo
echo "== fixture =="
python3 "$HERE/make_fixture.py" "$T/fixture.cif" >/dev/null
printf 'structure_id\tpdb_id\tchain\tcif_path\trevert\tassembly\nfake_all\tFAKE\tall\t%s\tauto\t\n' \
    "$T/fixture.cif" > "$T/structures.tsv"
ok "fixture written"

echo
echo "== find_gaps =="
python3 "$RG/bin/find_gaps.py" --samplesheet "$T/structures.tsv" \
    --out "$T/gaps.tsv" > "$T/find.log"
is "4 runs found" "$(( $(wc -l < "$T/gaps.tsv") - 1 ))" 4
is "2 eligible"   "$(awk -F'\t' 'NR>1 && $15=="yes"' "$T/gaps.tsv" | wc -l | tr -d ' ')" 2
# Author numbering must come out of the deposited scheme, not from the SEQRES
# index: the fixture offsets them by 100 on purpose.
has "gap reported in author numbering (106-108)" $'106\t108' "$T/gaps.tsv"
has "missing sequence read from SEQRES (STV)" "STV" "$T/gaps.tsv"
is "DNA chain not reported" "$(awk -F'\t' 'NR>1 && $3=="B"' "$T/gaps.tsv" | wc -l | tr -d ' ')" 0
is "terminal runs marked ineligible" \
   "$(awk -F'\t' 'NR>1 && $14!="internal" && $15=="yes"' "$T/gaps.tsv" | wc -l | tr -d ' ')" 0

python3 - "$T" <<'PY'
import sys, gemmi
from pathlib import Path
T = Path(sys.argv[1])
d = gemmi.cif.read(str(T / "fixture.cif")); b = d.sole_block()
b.find_mmcif_category("_pdbx_poly_seq_scheme.").erase()
d.write_file(str(T / "noscheme.cif"))
d2 = gemmi.cif.read(str(T / "fixture.cif")); b2 = d2.sole_block()
t = b2.find("_pdbx_poly_seq_scheme.", ["pdb_seq_num", "pdb_mon_id", "auth_mon_id"])
for r in t:
    if r.str(0) == "106":
        r[1] = "SER"; r[2] = "SER"
d2.write_file(str(T / "liar.cif"))
PY
for name in noscheme liar; do
    printf 'structure_id\tpdb_id\tchain\tcif_path\n%s\tFAKE\tA\t%s\n' \
        "$name" "$T/$name.cif" > "$T/s_$name.tsv"
done
printf 'structure_id\tpdb_id\tchain\tcif_path\ndna\tFAKE\tB\t%s\n' \
    "$T/fixture.cif" > "$T/s_dna.tsv"
dies python3 "$RG/bin/find_gaps.py" --samplesheet "$T/s_noscheme.tsv" --out "$T/x.tsv"
dies python3 "$RG/bin/find_gaps.py" --samplesheet "$T/s_liar.tsv" --out "$T/x.tsv"
dies python3 "$RG/bin/find_gaps.py" --samplesheet "$T/s_dna.tsv" --out "$T/x.tsv"

echo
echo "== make_inputs =="
# --min-fixed 6 because the fixture chain is 18 residues long and its merged
# loop covers 12 of them. The production default of 8 is unchanged; this is the
# fixture being small, not the threshold being relaxed. The negative case below
# still trips the same guard.
python3 "$RG/bin/make_inputs.py" --gaps "$T/gaps.tsv" --out-dir "$T/jobs" \
    --rosetta-db /DB --loopmodel-bin loopmodel --remodel-bin remodel \
    --nstruct 50 --chunk 25 --min-fixed 6 > "$T/make.log"
W="$T/jobs/fake_all"
is "18 full-pose residues" "$(( $(wc -l < "$W/numbering.tsv") - 1 ))" 18
is "5 residues marked modelled" \
   "$(awk -F'\t' 'NR>1 && $7==1' "$W/numbering.tsv" | wc -l | tr -d ' ')" 5
is "expected sequence" "$(cat "$W/expected_sequence.txt")" "GKLMSTVWYACDEFGHIK"
is "5 insertion lines in the blueprint" \
   "$(grep -c '^0 x ' "$W/insert.blueprint")" 5
is "blueprint covers every observed residue" \
   "$(grep -vc '^0 x ' "$W/insert.blueprint")" 13
# The fixture's two gaps are close enough that their flanked spans overlap, so
# they must become ONE loop rather than two overlapping LOOP lines.
is "overlapping spans merged into one loop" "$(wc -l < "$W/loops/all.loop" | tr -d ' ')" 1
has "loop in full-pose numbering" "LOOP 3 14 0 0 1" "$W/loops/all.loop"
has "NGK flag set present" "loops:kic_omega_sampling" "$W/flags.ngk"
has "score file requested" "out:file:scorefile" "$W/flags.ngk"
# A flank big enough to reach a modelled residue must be refused, not silently
# anchored on a position Rosetta invented; and a flank that swallows the whole
# chain must be refused too, because it leaves the drift check and the spread
# superposition with no fixed residues and they would then pass vacuously.
dies python3 "$RG/bin/make_inputs.py" --gaps "$T/gaps.tsv" \
    --out-dir "$T/jobs_bigflank" --rosetta-db /DB --loopmodel-bin loopmodel \
    --remodel-bin remodel --flank 0
dies python3 "$RG/bin/make_inputs.py" --gaps "$T/gaps.tsv" \
    --out-dir "$T/jobs_toolong" --rosetta-db /DB --loopmodel-bin loopmodel \
    --remodel-bin remodel --flank 40 --min-fixed 6

echo
echo "== check_models, stage insert =="
for defect in none break sequence drift; do
    rm -f "${W:?}"/[0-9]*.pdb "${W:?}"/insert_selected.pdb
    python3 "$HERE/make_fake_models.py" --work-dir "$W" --mode insert -n 3 \
        --defect "$defect" >/dev/null
    python3 "$RG/bin/check_models.py" --work-dir "$W" --stage insert \
        >"$T/insert_$defect.log" 2>&1 || true
    verdict=$(awk -F'\t' 'NR==2{print $2}' "$W/insert_report.tsv")
    if [ "$defect" = none ]; then
        is "clean candidate accepted" "$verdict" "yes"
    else
        is "$defect candidate rejected" "$verdict" "no"
    fi
done
has "break reported as a C-N distance" "C-N" "$W/insert_report.tsv" || true
rm -f "${W:?}"/[0-9]*.pdb
python3 "$HERE/make_fake_models.py" --work-dir "$W" --mode insert -n 1 >/dev/null
python3 "$RG/bin/check_models.py" --work-dir "$W" --stage insert >/dev/null 2>&1
[ -f "$W/insert_selected.pdb" ] && ok "insert_selected.pdb written" \
    || bad "insert_selected.pdb missing"

echo
echo "== check_models, stage ngk =="
rm -rf "${W:?}"/models
python3 "$HERE/make_fake_models.py" --work-dir "$W" --mode ngk -n 6 \
    --amplitude 0.0 >/dev/null
python3 "$RG/bin/check_models.py" --work-dir "$W" --stage ngk \
    > "$T/ngk.log" 2>&1
is "6 models valid" "$(awk -F'\t' 'NR>1 && $2=="yes"' "$W/ngk_report.tsv" | wc -l | tr -d ' ')" 6
has "scores read from beside the models" "-106.000" "$W/ngk_report.tsv"
is "both gaps in the spread report" "$(( $(wc -l < "$W/spread_report.tsv") - 1 ))" 2
has "identical models report zero spread" "0.00" "$W/spread_report.tsv"
# An empty model set is a failure, not an empty result.
rm -rf "${W:?}/models" "${W:?}/runs"
dies python3 "$RG/bin/check_models.py" --work-dir "$W" --stage ngk
python3 "$HERE/make_fake_models.py" --work-dir "$W" --mode ngk -n 6 \
    --amplitude 0.0 >/dev/null
python3 "$RG/bin/check_models.py" --work-dir "$W" --stage ngk >/dev/null 2>&1

echo
echo "== spread arithmetic =="
if python3 "$HERE/test_spread.py" > "$T/spread.log" 2>&1; then
    ok "all analytic spread cases"
else
    bad "spread arithmetic"; tail -8 "$T/spread.log"
fi

echo
echo "== promote =="
python3 "$RG/bin/promote.py" --work-dir "$W" --model "$W/models/ngk_0006.pdb" \
    --out "$T/out/fake.repaired.cif" > "$T/promote.log" 2>&1
has "seq_dif preserved" "_struct_ref_seq_dif." "$T/out/fake.repaired.cif"
has "assembly preserved" "_pdbx_struct_assembly." "$T/out/fake.repaired.cif"
has "scheme updated" "marked 5 residue(s) observed" "$T/promote.log"
# An unchecked or failing model must not be promotable.
dies python3 "$RG/bin/promote.py" --work-dir "$W" \
    --model "$W/models/does_not_exist.pdb" --out "$T/out/nope.cif"

printf 'structure_id\tpdb_id\tchain\tcif_path\trepaired\tFAKE\tA\t%s\n' \
    "$T/out/fake.repaired.cif" >/dev/null
printf 'structure_id\tpdb_id\tchain\tcif_path\nrepaired\tFAKE\tA\t%s\n' \
    "$T/out/fake.repaired.cif" > "$T/s_repaired.tsv"
python3 "$RG/bin/find_gaps.py" --samplesheet "$T/s_repaired.tsv" \
    --out "$T/gaps_after.tsv" > "$T/after.log"
is "no internal gap left" \
   "$(awk -F'\t' 'NR>1 && $14=="internal"' "$T/gaps_after.tsv" | wc -l | tr -d ' ')" 0
is "two terminal runs remain" \
   "$(( $(wc -l < "$T/gaps_after.tsv") - 1 ))" 2

echo
echo "== prepare_structures accepts the promoted file =="
printf '>FAKE01 synthetic reference\nAGKMMSTVWYACDEFGHIKL\n' > "$T/reference.fasta"
if python3 "$REPO/bin/prepare_structures.py" --samplesheet "$T/s_repaired.tsv" \
        --reference "$T/reference.fasta" --out-dir "$T/prep" > "$T/prep.log" 2>&1; then
    has "18 of 20 residues observed" "obs=18/20" "$T/prep.log"
    has "one substitution reverted" "reverted=1" "$T/prep.log"
    crosscheck=$(python3 -c "
import json,sys
print(json.load(open('$T/prep/structure_preparation.json'))['structures'][0].get('crosscheck_agrees'))")
    is "deposited substitution record still cross-checks" "$crosscheck" "True"
else
    bad "prepare_structures.py rejected the promoted file"; tail -5 "$T/prep.log"
fi

echo
echo "-- $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ] || exit 1
