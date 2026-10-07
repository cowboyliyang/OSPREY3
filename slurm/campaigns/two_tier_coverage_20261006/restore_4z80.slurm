#!/bin/bash
#SBATCH --job-name=restore-4z80-p0
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/restore_4z80_p0_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/restore_4z80_p0_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/two_tier_coverage_20261006/restore_4z80.py
