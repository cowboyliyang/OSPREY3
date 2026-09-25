#!/bin/bash
#SBATCH --job-name=ps_pair38_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 "$PAIR38_ROOT/source/pair38_20260921.py" summarize --root "$PAIR38_ROOT"
