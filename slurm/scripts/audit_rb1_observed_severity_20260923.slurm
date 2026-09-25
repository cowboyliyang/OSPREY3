#!/bin/bash
#SBATCH --job-name=rb1_s0_audit
#SBATCH --account=grisman
#SBATCH --partition=compsci,grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/usr/xtmp/lz280
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/audit_rb1_observed_severity_20260923.py
