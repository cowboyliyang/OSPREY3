#!/bin/bash
#SBATCH --job-name=pair38_rb24_verify
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/verify_pair38_rb24_20260922.py
