#!/bin/bash
#SBATCH --job-name=rb24_final_s0_audit
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
# Output estimate: 4 audit tables/JSON/Markdown, comfortably below 5 MiB total.
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/audit_rb24_final_severity_20260924.py
