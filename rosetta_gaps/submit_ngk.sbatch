#!/bin/bash
# Stage 2: next-generation KIC resamples every inserted span.
# One array task per row of jobs.tsv. Submit from the project directory on
# BeeGFS, after stage 1 has produced insert_selected.pdb:
#
#   sbatch --array=1-$(( $(wc -l < rosetta_gaps/work/jobs.tsv) - 1 )) \
#          rosetta_gaps/submit_ngk.sbatch rosetta_gaps/work containers/rosetta.sif
#
#SBATCH --job-name=rosetta-ngk
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/ngk_%A_%a.out
#SBATCH --error=logs/ngk_%A_%a.err
#
# 12G because -loops:kic_rama2b is documented to raise memory use to about 6G,
# and KIC-with-fragments also holds a fragment library. No time limit, by
# project rule.
#
# UNVERIFIED: nothing in this repository has run Rosetta. Read MaxRSS off the
# first array with `sacct -j <jobid> --format=JobID,MaxRSS,State` and lower it.
set -euo pipefail

ROOT="${1:?usage: submit_ngk.sbatch <jobs root> <rosetta.sif>}"
SIF="${2:?usage: submit_ngk.sbatch <jobs root> <rosetta.sif>}"
TASK="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"

case "$PWD" in
    /home/*/isilon/*) echo "Launch from BeeGFS, not NFS." >&2; exit 1 ;;
esac

JOBS="$ROOT/jobs.tsv"
[ -f "$JOBS" ] || { echo "no $JOBS; run bin/make_inputs.py first" >&2; exit 1; }
[ -f "$SIF" ] || { echo "no image $SIF" >&2; exit 1; }

read -r WORK FLAGS STYLE NSTRUCT ROW_TASK < <(
    awk -F'\t' -v want="$TASK" '
        NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        NR - 1 == want { print $col["work_dir"], $col["flags"],
                               $col["style"], $col["nstruct"], $col["task"] }
    ' "$JOBS")
[ -n "${WORK:-}" ] || { echo "array task $TASK has no row in $JOBS" >&2; exit 1; }

PROJECT="$PWD"
module load singularity
mkdir -p logs

# Every task runs in its own directory. Rosetta's output names restart at
# _0001 in each task, so a shared directory would have tasks overwriting each
# other's models - and the score file would then describe the wrong structure.
RUN="$WORK/runs/t${TASK}_${STYLE}"
mkdir -p "$RUN"
cd "$RUN"
for item in insert_selected.pdb loops "$FLAGS"; do
    ln -sfn "$PROJECT/$WORK/$item" "./$(basename "$item")"
done

[ -e insert_selected.pdb ] || { echo "stage 1 produced no insert_selected.pdb \
in $WORK: run submit_insert.sbatch and check its insert_report.tsv" >&2; exit 1; }

LOOPMODEL="${LOOPMODEL_BIN:-loopmodel}"

# A seed that depends on the array index: without -run:constant_seed Rosetta
# seeds from the system device and the run is not reproducible, and with a
# constant seed every task would produce the SAME models and the ensemble
# would be an illusion.
SEED=$(( 1111111 + TASK ))

ARGS=( "@$(basename "$FLAGS")" -nstruct "$NSTRUCT"
       -run:constant_seed -run:jran "$SEED" )

if [ "$STYLE" = "kicfrag" ]; then
    # UNVERIFIED: KIC with fragments needs fragment files, and where they come
    # from depends on the image. Set FRAG_FILES and FRAG_SIZES before
    # submitting a kicfrag task, e.g.
    #   export FRAG_FILES=/path/frags.200.9mers FRAG_SIZES=9
    : "${FRAG_FILES:?kicfrag needs FRAG_FILES (see the comment above)}"
    : "${FRAG_SIZES:?kicfrag needs FRAG_SIZES (see the comment above)}"
    ARGS+=( -loops:frag_files "$FRAG_FILES" -loops:frag_sizes "$FRAG_SIZES" )
fi

echo "task $TASK: $STYLE nstruct=$NSTRUCT seed=$SEED in $RUN"
singularity exec --cleanenv -B "$PROJECT" "$SIF" "$LOOPMODEL" "${ARGS[@]}"

echo "task $TASK done. When the whole array has finished, run:"
echo "  python3 rosetta_gaps/bin/check_models.py --work-dir $WORK --stage ngk"
