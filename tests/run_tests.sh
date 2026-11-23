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
