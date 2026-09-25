#!/bin/bash
#SBATCH --job-name=final_evalue_audit
#SBATCH --account=grisman
#SBATCH --partition=compsci,grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/usr/xtmp/lz280
# Five small audit tables/reports, estimated below 15 MiB. Frozen inputs are read only.
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/audit_final_evalue_20260924.py
