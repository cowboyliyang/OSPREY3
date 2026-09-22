#!/bin/bash
#SBATCH --job-name=h200_active162
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${PREP_ROOT:?Validated 25-design prep}"
: "${MARKSTAR_ARRAY:?Submitted 25-task array}"
ROOT=/usr/xtmp/lz280/markstar_backfill_20260922/handoff_${SLURM_JOB_ID}
mkdir -p "$ROOT/tmp" "$ROOT/scripts"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/handoff_active162_20260922.py "$ROOT/scripts/handoff.py"
cp "${BASH_SOURCE[0]}" "$ROOT/scripts/handoff.slurm"
python3 "$ROOT/scripts/handoff.py" "$PREP_ROOT" "$ROOT" "$MARKSTAR_ARRAY"
