#!/bin/bash
#SBATCH --job-name=ps-auto-summary
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export AUTO_BENCH_SUMMARY=/usr/xtmp/lz280/packstar_auto_preflight/summary_${SLURM_JOB_ID}
export PYTHONPYCACHEPREFIX="$AUTO_BENCH_SUMMARY/pycache"
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_auto_preflight_benefit.py "$@"
