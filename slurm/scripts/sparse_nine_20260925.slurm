#!/bin/bash
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --job-name=sp9_prepare
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/packstar_sparse_nine_20260925/launch_$SLURM_JOB_ID
mkdir -p "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/sparse_nine_20260925.py "$ROOT"
