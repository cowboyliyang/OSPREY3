#!/bin/bash
#SBATCH --job-name=frontier-extension7-prep
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/frontier_extension7_prep_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/frontier_extension7_prep_%j.err
set -euo pipefail
export TMPDIR=/usr/xtmp/lz280/frontier_extension7_20261006/tmp_prep_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/frontier_extension7_20261006/prepare.py
