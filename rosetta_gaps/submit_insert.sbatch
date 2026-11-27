#!/bin/bash
# Stage 1: RosettaRemodel inserts the missing residues and closes the chain.
# One array task per structure. Submit from the project directory on BeeGFS:
#
#   sbatch --array=1-$(( $(awk -F'\t' 'NR>1{print $6}' rosetta_gaps/work/jobs.tsv \
#            | sort -u | wc -l) )) \
#          rosetta_gaps/submit_insert.sbatch rosetta_gaps/work containers/rosetta.sif
#
#SBATCH --job-name=rosetta-insert
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/insert_%A_%a.out
#SBATCH --error=logs/insert_%A_%a.err
#
# Rosetta is single-process here, so one task and one CPU. No time limit is
# set, by project rule.
#
# UNVERIFIED: 8G is a guess. Nothing in this repository has run Rosetta.
# Read the real figure off the first run with
#   sacct -j <jobid> --format=JobID,MaxRSS,State
# and lower it.
set -euo pipefail

ROOT="${1:?usage: submit_insert.sbatch <jobs root> <rosetta.sif>}"
SIF="${2:?usage: submit_insert.sbatch <jobs root> <rosetta.sif>}"
TASK="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"

case "$PWD" in
    /home/*/isilon/*) echo "Launch from BeeGFS, not NFS." >&2; exit 1 ;;
esac

JOBS="$ROOT/jobs.tsv"
[ -f "$JOBS" ] || { echo "no $JOBS; run bin/make_inputs.py first" >&2; exit 1; }
[ -f "$SIF" ] || { echo "no image $SIF" >&2; exit 1; }

# Unique work directories, in a stable order, so a task id always maps to the
# same structure across resubmissions.
WORK=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="work_dir") c=i; next}
                   {print $c}' "$JOBS" | sort -u | sed -n "${TASK}p")
[ -n "$WORK" ] || { echo "array task $TASK has no work_dir in $JOBS" >&2; exit 1; }

PROJECT="$PWD"
module load singularity
mkdir -p logs

echo "task $TASK: $WORK"
cd "$WORK"

# The binary name comes from bin/probe_container.sh, not from a guess.
REMODEL="${REMODEL_BIN:-remodel}"

singularity exec --cleanenv -B "$PROJECT" "$SIF" \
    "$REMODEL" @flags.remodel \
    -run:constant_seed -run:jran $(( 1111111 + TASK ))

# Remodel's job-distributor output (<input>_0001.pdb) is not the model to use
# when num_trajectory > 1; the trajectory outputs are 1.pdb, 2.pdb, ... This
# check reads those and writes insert_selected.pdb.
python3 "$PROJECT/rosetta_gaps/bin/check_models.py" \
    --work-dir . --stage insert

echo "task $TASK done"
