#!/bin/bash
#SBATCH --job-name=packstar_forward_summary
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${COMPARISON_JOB:?comparison array id required}"
export TMPDIR=/usr/xtmp/lz280/packstar_forward_comparison_20260913/A${COMPARISON_JOB}/tmp-summary
mkdir -p "$TMPDIR"
python3 "${BUILD_ROOT:?tested build root required}/source/slurm/scripts/summarize_forward_triples_20260913.py" \
    /usr/xtmp/lz280/packstar_forward_comparison_20260913/A${COMPARISON_JOB}
