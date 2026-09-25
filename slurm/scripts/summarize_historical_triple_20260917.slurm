#!/bin/bash
#SBATCH --job-name=ps_old_triple_summary
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
: "${OLD_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_historical_triple_20260917/A${OLD_JOB}
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_historical_triple_20260917.py "$ROOT/summarize.py"
python3 "$ROOT/summarize.py" "$OLD_JOB"
