#!/bin/bash
#SBATCH --job-name=ps_million_merge
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
AUDIT="${RESULT_ROOT:?}/merge_script_$SLURM_JOB_ID"
mkdir -p "$AUDIT"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/merge_million_census_20260918.py "$AUDIT/"
export PYTHONPYCACHEPREFIX="$AUDIT/pycache"
python3 "$AUDIT/merge_million_census_20260918.py" "$RESULT_ROOT" "${ARRAY_ID:?}" \
 "${BUILD_ROOT:?}/screen/${SYSTEM:?}.tsv" "${BLOCK:?}"
