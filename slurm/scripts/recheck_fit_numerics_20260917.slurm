#!/bin/bash
#SBATCH --job-name=ps_fit_recheck
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/recheck_${SLURM_JOB_ID}
mkdir -p "$ROOT"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/check_proposal_fit_numerics_20260917.py "$ROOT/check.py"
cp "${BASH_SOURCE[0]}" "$ROOT/recheck.slurm"
export PYTHONPYCACHEPREFIX="$ROOT/pycache"
for CASE in 1gwc 2xgy; do
    python3 "$ROOT/check.py" \
        /usr/xtmp/lz280/packstar_fit_speed_20260917/refine_12626643/threads16 \
        /usr/xtmp/lz280/packstar_fit_speed_20260917/kernel_12626677/shifted \
        "$CASE" "$ROOT/${CASE}_check.json"
done
touch "$ROOT/PASSED"
