#!/bin/bash
#SBATCH --job-name=restart-packstar-gpu
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/restart_packstar_gpu_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/restart_packstar_gpu_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/gpu_memory_retry_20261006/restart.py
