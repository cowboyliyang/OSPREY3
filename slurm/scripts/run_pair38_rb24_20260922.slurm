#!/bin/bash
#SBATCH --job-name=pair38_rb24
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -f "${PAIR38_ROOT:?}/READY"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
python3 "$PAIR38_ROOT/source/pair38_rb24_20260922.py" run --root "$PAIR38_ROOT" --index "${CASE_INDEX:?}"
