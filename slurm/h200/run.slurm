#!/bin/bash
# Supply account, partition, CPUs, memory, GPUs, array and logs to sbatch.
#SBATCH --job-name=packstar-h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${BUILD_ROOT:?Validated build directory}"
test -f "$BUILD_ROOT/READY"
export PYTHONPYCACHEPREFIX="$BUILD_ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 "$BUILD_ROOT/source/slurm/h200/run.py" "$@"
