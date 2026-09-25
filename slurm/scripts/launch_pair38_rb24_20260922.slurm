#!/bin/bash
#SBATCH --job-name=pair38_rb24_launch
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/pair38_rb24_20260922.py launch
