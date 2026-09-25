#!/bin/bash
#SBATCH --job-name=markstar_flex_add5_128g
#SBATCH --account=grisman
#SBATCH --partition=compsci,grisman
#SBATCH --exclude=fennario-[01-06]
#SBATCH --array=24,25
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
PREP_ROOT=/usr/xtmp/lz280/packstar_flex_frontier12_20260918/prep_12632172
REPO=/home/users/lz280/OSPREY3-fresh-packstar
test -s "$PREP_ROOT/READY"
TASK_ROOT=/usr/xtmp/lz280/packstar_flex_frontier12_20260918/A${SLURM_ARRAY_JOB_ID}/T${SLURM_ARRAY_TASK_ID}
test ! -e "$TASK_ROOT/run_manifest.json"
mkdir -p "$TASK_ROOT/tmp"
export TMPDIR="$TASK_ROOT/tmp" PYTHONPYCACHEPREFIX="$TASK_ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
cp "${BASH_SOURCE[0]}" "$TASK_ROOT/submitted.slurm"
cp "$REPO/slurm/scripts/execute_markstar_flex_add5_128g_20260918.py" "$TASK_ROOT/execute.py"
python3 "$TASK_ROOT/execute.py" "$PREP_ROOT" "$TASK_ROOT" "$SLURM_ARRAY_TASK_ID"
