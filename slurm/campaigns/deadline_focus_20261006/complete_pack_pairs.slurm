#!/bin/bash
#SBATCH --job-name=complete-pack-pairs
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/complete_pack_pairs_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/complete_pack_pairs_%j.err
set -euo pipefail
export TMPDIR=/usr/xtmp/lz280/markstar_deadline_focus_20261006/tmp_pairs_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/deadline_focus_20261006/complete_pack_pairs.py
