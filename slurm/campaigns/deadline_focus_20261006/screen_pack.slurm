#!/bin/bash
#SBATCH --account=grisman
#SBATCH --job-name=deadline-pack-screen
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/deadline_pack_screen_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/deadline_pack_screen_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 "$EXTENSION_PREP/screen_pack.py"
