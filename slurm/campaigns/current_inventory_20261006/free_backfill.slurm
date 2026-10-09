#!/bin/bash
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --job-name=mark-backfill-controller
#SBATCH --cpus-per-task=1
#SBATCH --mem=128M
#SBATCH --time=90-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/mark_backfill_controller_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/mark_backfill_controller_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
exec python3 "${MARKSTAR_BACKFILL_SCRIPT:-/home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/free_backfill.py}"
