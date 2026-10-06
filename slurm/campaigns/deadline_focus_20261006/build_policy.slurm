#!/bin/bash
#SBATCH --account=grisman
#SBATCH --job-name=deadline-focus-registry
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/deadline_focus_registry_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/deadline_focus_registry_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/deadline_focus_20261006/build_policy.py
