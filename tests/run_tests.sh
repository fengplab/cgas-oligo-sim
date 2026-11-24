#!/bin/bash
# End-to-end test of the visualisation branch on synthetic fixtures.
# Runs the bin/ scripts directly: it does NOT exercise Nextflow, Slurm or
# Singularity. Those are cluster-side and untested from here.
set -euo pipefail
cd "$(dirname "$0")/.."

rm -rf tests/fixtures tests/out
mkdir -p tests/fixtures tests/out/figs tests/out/report/figures

python3 bin/make_synthetic.py --outdir tests/fixtures

python3 bin/call_sites.py \
    --occupancy tests/fixtures/occupancy.tsv \
    --null-occupancy tests/fixtures/null_occupancy.tsv \
    --null-percentile 90 --dwell-min 0.05 --bootstrap-n 2000 --seed 20261120 \
    --replicate-min "" \
    --out-sites tests/out/sites.tsv --out-profile tests/out/profile.tsv \
    --out-tiers tests/out/tiers.tsv --out-summary tests/out/summary.json

# Fixture-only: rekey synthetic poses onto the called site ids, matched by
echo
echo "--- environment.yml matches the code's imports ---"
python3 tests/check_environment.py

echo
echo "--- prepare_structures: find crystallisation mutations by ALIGNMENT ---"
# The fixture carries K186N and L195R relative to its reference. No mutation
# table is given to the script; it must recover both from the alignment, and
# must refuse a reference that is the wrong protein.
python3 tests/make_structure_fixture.py
python3 bin/prepare_structures.py --samplesheet tests/fixtures/structures.tsv \
    --reference tests/fixtures/structure_reference.fasta --out-dir tests/out/struct
python3 - <<'PYCHK'
import sys
import gemmi
import pandas as pd
rv = pd.read_csv("tests/out/struct/wt_reversions.tsv", sep="\t")
got = {(int(r.auth_seq), r["from"], r["to"]) for _, r in rv.iterrows()}
want = {(186, "ASN(N)", "LYS(K)"), (195, "ARG(R)", "LEU(L)")}
if got != want:
    sys.exit(f"FAIL: reversions {sorted(got)} != {sorted(want)}")
st = gemmi.read_structure("tests/out/struct/structures/FIXWT.wt.cif")
for res in st[0]["A"]:
    if res.seqid.num in (186, 195):
        names = sorted(a.name for a in res)
        if names != ["C", "CA", "CB", "N", "O"]:
            sys.exit(f"FAIL: residue {res.seqid.num} kept {names}; strip policy "
                     "must leave only backbone + CB so no contact is derived "
                     "from a sidechain the wild-type residue lacks")
print("  OK K186->LYS, L195->LEU, sidechains stripped to backbone+CB")
PYCHK

echo "--- prepare_structures: the wrong reference must FAIL, not revert blindly ---"
printf '>unrelated\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG\n' \
    > tests/out/wrong_reference.fasta
if python3 bin/prepare_structures.py --samplesheet tests/fixtures/structures.tsv \
    --reference tests/out/wrong_reference.fasta --out-dir tests/out/struct_wrong 2>/dev/null; then
    echo "FAIL: an unrelated reference sequence was accepted"; exit 1
fi
echo "  OK rejected"

echo
echo "--- enumerate_arms: PWMcos arms appear ONLY when motifs are passed ---"
printf 'name\nrank1\nrank2\nrank3\n' > tests/out/arm_panel.tsv
printf 'structure_id\nS1\nS2\n' > tests/out/arm_structures2.tsv
printf 'structure_id\nS1\n' > tests/out/arm_structures1.tsv
printf 'family\nfamA\nfamB\n' > tests/out/arm_motifs.tsv

python3 bin/enumerate_arms.py --panel tests/out/arm_panel.tsv \
    --structures tests/out/arm_structures1.tsv --partners TFAM \
    --replicates 10 --shuffles 4 --null-replicates 3 \
    --out-jobs tests/out/jobs_nomotif.tsv --out-arms tests/out/arms_nomotif.tsv \
    --out-summary tests/out/arms_nomotif.json

python3 bin/enumerate_arms.py --panel tests/out/arm_panel.tsv \
    --structures tests/out/arm_structures2.tsv --motifs tests/out/arm_motifs.tsv \
    --partners TFAM,HMGB1 --replicates 10 --shuffles 4 --null-replicates 3 \
    --out-jobs tests/out/jobs_full.tsv --out-arms tests/out/arms_full.tsv \
    --out-summary tests/out/arms_full.json
python3 - <<'PYCHK'
import json, sys
import pandas as pd
a = json.load(open("tests/out/arms_nomotif.json"))
if a["pwmcos_arms"] != 0 or not a["pwmcos_skipped_no_motifs"]:
    sys.exit("FAIL: PWMcos arms generated with no motifs passed")
if a["arms"] != ["baseline", "partner:TFAM"]:
    sys.exit(f"FAIL: unexpected arms without motifs: {a['arms']}")
b = json.load(open("tests/out/arms_full.json"))
if b["pwmcos_arms"] != 2 or not b["multibasin_arm"]:
    sys.exit(f"FAIL: expected 2 PWMcos arms and a multibasin arm, got {b}")
j = pd.read_csv("tests/out/jobs_full.tsv", sep="\t")
if j.job_id.duplicated().any():
    sys.exit("FAIL: duplicate job_id would make jobs overwrite each other")
if j.seed.nunique() != len(j):
    sys.exit("FAIL: seeds are not unique per job")
print(f"  OK {a['n_arms']} arms without motifs, {b['n_arms']} with; "
      f"{len(j)} unique jobs")
PYCHK

echo
echo "--- extract_occupancy: synthetic trajectory with a known register ---"
# The fixture puts one cGAS copy over nucleotides 28-45, so the register the
# extractor reports must be their midpoint, 36. A 10 A cutoff in 3D also
# catches beads ~2.5 nt beyond each end, which widens the CONTACT window but
# must not move the register.
python3 tests/make_occupancy_fixture.py --out-prefix tests/out/occfix
python3 bin/extract_occupancy.py \
    --trajectory tests/out/occfix.npz --bead-map tests/out/occfix.bead_map.tsv \
    --sequence-id FIX_A --replicate 1 --arm baseline \
    --out-occupancy tests/out/occfix.occupancy.tsv \
    --out-contacts tests/out/occfix.contacts.tsv \
    --out-orientation tests/out/occfix.orientation.tsv \
    --out-summary tests/out/occfix.summary.json
python3 - <<'PYCHK'
import json, sys
import pandas as pd
occ = pd.read_csv("tests/out/occfix.occupancy.tsv", sep="\t")
top = occ.loc[occ.dwell_frac.idxmax()]
if int(top.position) != 36:
    sys.exit(f"FAIL: register {int(top.position)} != 36 (window 28-45 midpoint)")
if not 0.45 < top.dwell_frac < 0.75:
    sys.exit(f"FAIL: dwell_frac {top.dwell_frac:.3f} outside 0.45-0.75 for a 0.6 fixture")
s = json.load(open("tests/out/occfix.summary.json"))
if s["frames_with_binding_frac"] > 1.0:
    sys.exit("FAIL: frames_with_binding_frac must never exceed 1")
print(f"  OK register=36 dwell={top.dwell_frac:.3f}")
PYCHK

echo "--- extract_occupancy: a replicate that never binds is a VALID result ---"
python3 tests/make_occupancy_fixture.py --out-prefix tests/out/occnb --bound-frac 0.0
python3 bin/extract_occupancy.py \
    --trajectory tests/out/occnb.npz --bead-map tests/out/occnb.bead_map.tsv \
    --sequence-id FIX_NB --replicate 1 \
    --out-occupancy tests/out/occnb.occupancy.tsv \
    --out-contacts tests/out/occnb.contacts.tsv \
    --out-orientation tests/out/occnb.orientation.tsv \
    --out-summary tests/out/occnb.summary.json
python3 -c "
import json,sys,pandas as pd
assert json.load(open('tests/out/occnb.summary.json'))['no_binding_events'] is True
assert pd.read_csv('tests/out/occnb.occupancy.tsv',sep='\t').dwell_frac.max() == 0.0
print('  OK zero occupancy, exit 0')"

echo "--- extract_occupancy: bead-map mismatch must FAIL ---"
head -250 tests/out/occfix.bead_map.tsv > tests/out/occfix.short.tsv
if python3 bin/extract_occupancy.py \
    --trajectory tests/out/occfix.npz --bead-map tests/out/occfix.short.tsv \
    --sequence-id X --replicate 1 --out-occupancy /dev/null --out-contacts /dev/null \
    --out-orientation /dev/null --out-summary /dev/null 2>/dev/null; then
    echo "FAIL: a bead map shorter than the trajectory was accepted"; exit 1
fi
echo "  OK rejected"

# position. The real pipeline has no such step - pose extraction runs
# downstream of CALL_SITES and is keyed on site_id directly.
python3 tests/rekey_fixture_poses.py

while IFS=$'\t' read -r seq sid npz; do
    [ "$seq" = "sequence_id" ] && continue
    python3 bin/cluster_poses.py --npz "$npz" --site-id "$sid" --sequence-id "$seq" \
        --contacts tests/fixtures/contacts.tsv --rmsd-cutoff 3.5 --min-frames 50 \
        --out-pdb "tests/out/${sid}.pdb" \
        --out-clusters "tests/out/${sid}_clusters.tsv" \
        --out-summary "tests/out/${sid}_cluster.json"
done < tests/fixtures/poses.tsv

python3 bin/make_figures.py --profile tests/out/profile.tsv --sites tests/out/sites.tsv \
    --occupancy tests/fixtures/occupancy.tsv --contacts tests/fixtures/contacts.tsv \
    --outdir tests/out/figs

cp tests/out/figs/*.png tests/out/report/figures/
python3 bin/build_index.py --sites tests/out/sites.tsv \
    --cluster-summaries tests/out/*_cluster.json --pdbs tests/out/*.pdb \
    --figdir tests/out/figs --out tests/out/report/index.html \
    --run-params '{"null_percentile":90,"dwell_min":0.05,"replicate_min":"unset","rmsd_cutoff":3.5,"seed":20261120}'

echo
echo "PASS: full chain completed. Report at tests/out/report/index.html"
