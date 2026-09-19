#!/bin/bash
# New 25-design cohort; reuse the validated H200 build without editing it.
#SBATCH --job-name=packstar-h200-add4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${BUILD_ROOT:?Existing validated H200 build}"
: "${REPO:?Updated source checkout containing the add-on manifest and runner}"
: "${RESULT_ROOT:?Separate scratch result root for the four additions}"
test -f "$BUILD_ROOT/READY"
DRIVER_ROOT="$RESULT_ROOT/drivers/J${SLURM_JOB_ID}_T${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$RESULT_ROOT/drivers"
mkdir "$DRIVER_ROOT"
cp "$REPO/slurm/h200/run.py" "$DRIVER_ROOT/run.py"
cp "$REPO/slurm/h200/frontier_add4.tsv" "$DRIVER_ROOT/frontier_add4.tsv"
cp "${BASH_SOURCE[0]}" "$DRIVER_ROOT/submitted.slurm"
export PYTHONPYCACHEPREFIX="$DRIVER_ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 "$DRIVER_ROOT/run.py" --frontier-manifest "$DRIVER_ROOT/frontier_add4.tsv" "$@"
