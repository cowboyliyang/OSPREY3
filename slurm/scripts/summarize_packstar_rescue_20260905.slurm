#!/bin/bash
# Aggregate completed rescue-array artifacts.  Analysis is kept on Slurm per
# the workspace execution policy.
#SBATCH --job-name=ps_rescue_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
SELECTION_ROOT=${SELECTION_ROOT:?SELECTION_ROOT is required}
RESULT_ROOT=${RESULT_ROOT:?RESULT_ROOT is required}
OUT_ROOT=${OUT_ROOT:-$RESULT_ROOT/summary_${SLURM_JOB_ID}}
mkdir -p "$OUT_ROOT"
python3 "$REPO/slurm/scripts/summarize_packstar_rescue.py" \
    --candidates "$SELECTION_ROOT/candidates.tsv" \
    --result-root "$RESULT_ROOT" \
    --output "$OUT_ROOT/results.tsv"
