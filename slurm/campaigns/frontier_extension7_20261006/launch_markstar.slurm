#!/bin/bash
#SBATCH --job-name=markstar-extension-launch
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_extension_launch_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_extension_launch_%j.err
set -euo pipefail
: "${MARKSTAR_EXTENSION_PREP:?}"
export PYTHONDONTWRITEBYTECODE=1
python3 "$MARKSTAR_EXTENSION_PREP/launch_markstar.py"
