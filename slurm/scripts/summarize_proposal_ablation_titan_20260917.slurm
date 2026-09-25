#!/bin/bash
#SBATCH --job-name=ps_ablation_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${RESULT_ROOT:?}" "${BUILD_ROOT:?}"
mkdir -p "$RESULT_ROOT/tmp"
export TMPDIR="$RESULT_ROOT/tmp"
export PYTHONPYCACHEPREFIX="$RESULT_ROOT/pycache"
python3 "$BUILD_ROOT/source/slurm/scripts/summarize_proposal_ablation_titan_20260917.py" "$RESULT_ROOT"
