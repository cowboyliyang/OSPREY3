#!/bin/bash
#SBATCH --job-name=markstar-extension24-plan
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_extension24_plan_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_extension24_plan_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/frontier_extension7_20261006/plan_markstar.py
