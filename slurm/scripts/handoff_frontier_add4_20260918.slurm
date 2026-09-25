#!/bin/bash
#SBATCH --job-name=ms_add4_handoff
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${PREP_ROOT:?}"
: "${SMALL_JOB:?}"
: "${LARGE_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_flex_frontier16_20260918/handoff_${SLURM_JOB_ID}
export TMPDIR=/usr/xtmp/lz280/packstar_flex_frontier16_20260918/handoff_tmp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
export PYTHONPYCACHEPREFIX="$TMPDIR/pycache"
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/handoff_frontier_add4_20260918.py "$PREP_ROOT" "$ROOT" "$SMALL_JOB" "$LARGE_JOB"
