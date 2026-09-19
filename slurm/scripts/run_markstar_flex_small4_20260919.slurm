#!/bin/bash
#SBATCH --job-name=markstar_flex_small4
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --exclude=fennario-[01-06]
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
: "${PREP_ROOT:?}"
test -s "$PREP_ROOT/READY"
TASK_ROOT=/usr/xtmp/lz280/packstar_flex_frontier20_20260919/A${SLURM_ARRAY_JOB_ID}/T${SLURM_ARRAY_TASK_ID}
test ! -e "$TASK_ROOT/run_manifest.json"
mkdir -p "$TASK_ROOT/tmp"
export TMPDIR="$TASK_ROOT/tmp" PYTHONPYCACHEPREFIX="$TASK_ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
cp "${BASH_SOURCE[0]}" "$TASK_ROOT/submitted.slurm"
python3 "$PREP_ROOT/scripts/execute.py" "$PREP_ROOT" "$TASK_ROOT" "$SLURM_ARRAY_TASK_ID"
