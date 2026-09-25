#!/bin/bash
#SBATCH --job-name=ps_hardcase_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
SUMMARY_CODE=${SUMMARY_CODE:-$BUILD_ROOT/source/slurm/scripts/hardcase_ablation_20260921.py}
cp "$SUMMARY_CODE" "$RESULT_ROOT/summarizer_${SLURM_JOB_ID}.py"
sha256sum "$RESULT_ROOT/summarizer_${SLURM_JOB_ID}.py" > "$RESULT_ROOT/summarizer_${SLURM_JOB_ID}.sha256"
python3 "$RESULT_ROOT/summarizer_${SLURM_JOB_ID}.py" summarize \
    --build "$BUILD_ROOT" --root "$RESULT_ROOT"
