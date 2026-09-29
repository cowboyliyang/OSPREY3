#!/bin/bash
#SBATCH --job-name=ps-log10-postcheck
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_postcheck_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_postcheck_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${POSTCHECK_PACKAGE:?Frozen follow-up scripts}"
cd "$POSTCHECK_PACKAGE"
sha256sum -c SHA256SUMS
export FOLLOWUP_ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/postcheck_${SLURM_JOB_ID}
mkdir "$FOLLOWUP_ROOT"
mkdir "$FOLLOWUP_ROOT/tmp"
export TMPDIR="$FOLLOWUP_ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
python3 "$POSTCHECK_PACKAGE/log10_postcheck_20260927.py"
