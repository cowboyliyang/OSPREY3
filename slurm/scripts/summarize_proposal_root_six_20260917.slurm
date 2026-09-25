#!/bin/bash
#SBATCH --job-name=ps_proposal_root_summary
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${RUN_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_proposal_root_20260917/A${RUN_JOB}
mkdir -p "$ROOT/summary_tmp"
export TMPDIR="$ROOT/summary_tmp" PYTHONPYCACHEPREFIX="$ROOT/summary_pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/summarize.slurm"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_proposal_root_six_20260917.py "$ROOT/summarize.py"
python3 "$ROOT/summarize.py" "$RUN_JOB"
