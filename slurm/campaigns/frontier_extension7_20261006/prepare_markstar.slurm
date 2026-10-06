#!/bin/bash
#SBATCH --job-name=markstar-extension-prep
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_extension_prep_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_extension_prep_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/frontier_extension7_20261006/prepare_markstar.py
