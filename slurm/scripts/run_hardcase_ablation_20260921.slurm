#!/bin/bash
#SBATCH --job-name=ps_hardcase_ablation
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --exclusive
#SBATCH --gres=gpu:a5000:4
#SBATCH --array=0-11%6
#SBATCH --time=16:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:?}"
test -f "${BUILD_ROOT:?}/READY"
test -n "${RESULT_ROOT:?}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
python3 "$BUILD_ROOT/source/slurm/scripts/hardcase_ablation_20260921.py" run \
    --build "$BUILD_ROOT" --root "$RESULT_ROOT" --task "${SLURM_ARRAY_TASK_ID:-0}"
