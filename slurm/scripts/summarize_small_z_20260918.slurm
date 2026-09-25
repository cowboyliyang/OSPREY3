#!/bin/bash
#SBATCH --job-name=ps_z_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONPYCACHEPREFIX="${RESULT_ROOT:?}/pycache"
AUDIT="$RESULT_ROOT/analysis_$SLURM_JOB_ID"
mkdir -p "$AUDIT"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_small_z_20260918.py "$AUDIT/"
sha256sum "$AUDIT/summarize_small_z_20260918.py" > "$AUDIT/source.sha256"
python3 "$AUDIT/summarize_small_z_20260918.py" "$RESULT_ROOT" "${SUMMARY_MODE:-summary}"
