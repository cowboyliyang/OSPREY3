#!/bin/bash
#SBATCH --job-name=verify-markstar-continuation
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/verify_markstar_continuation_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/verify_markstar_continuation_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/deadline_focus_20261006/verify_cpu_continuation.py
