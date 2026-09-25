#!/bin/bash
#SBATCH --job-name=pair38_interim
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
# Small reports only, expected <5 MiB / 10 files; run data are read, not modified.
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/status_pair38_20260921.py
