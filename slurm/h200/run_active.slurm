#!/bin/bash
# Current 162-design manifest; reuse the validated H200 build.
#SBATCH --job-name=packstar-h200-active
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${BUILD_ROOT:?Existing validated H200 build}"
: "${REPO:?Updated checkout with frontier_active.tsv}"
: "${INPUT_ROOT:?Transferred active package}"
: "${RESULT_ROOT:?Scratch output root}"
test -f "$BUILD_ROOT/READY"
DRIVER_ROOT="$RESULT_ROOT/drivers/J${SLURM_JOB_ID}_T${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$RESULT_ROOT/drivers"
mkdir "$DRIVER_ROOT"
cp "$REPO/slurm/h200/run.py" "$DRIVER_ROOT/run.py"
cp "$REPO/slurm/h200/frontier_active.tsv" "$DRIVER_ROOT/frontier_active.tsv"
cp "${BASH_SOURCE[0]}" "$DRIVER_ROOT/submitted.slurm"
export PYTHONPYCACHEPREFIX="$DRIVER_ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 "$DRIVER_ROOT/run.py" --frontier-manifest "$DRIVER_ROOT/frontier_active.tsv" "$@"
