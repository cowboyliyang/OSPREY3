#!/bin/bash
#SBATCH --job-name=ps_mcov_summary
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
AUDIT="${RESULT_ROOT:?}/analysis_$SLURM_JOB_ID"
mkdir -p "$AUDIT"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_million_coverage_20260918.py "$AUDIT/"
export PYTHONPYCACHEPREFIX="$AUDIT/pycache"
python3 "$AUDIT/summarize_million_coverage_20260918.py" "$RESULT_ROOT" \
 "${REFERENCE_ROOT:-/usr/xtmp/lz280/packstar_million_census_20260918/run_12633933/merged_A12633943}"
