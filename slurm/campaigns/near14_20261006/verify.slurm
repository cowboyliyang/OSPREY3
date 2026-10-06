#!/bin/bash
#SBATCH --job-name=near14-verify
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/near14_verify_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/near14_verify_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/near14_20261006/verify.py
