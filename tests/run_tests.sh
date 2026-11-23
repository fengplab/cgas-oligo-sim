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
