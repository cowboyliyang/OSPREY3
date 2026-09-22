#!/bin/bash
#SBATCH --job-name=ms25_prepare
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/markstar_backfill_20260922/prep_${SLURM_JOB_ID}
REPO=/home/users/lz280/OSPREY3-fresh-packstar
mkdir -p "$ROOT/tmp" "$ROOT/scripts"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
cp "$REPO/slurm/scripts/prepare_backfill25_20260922.py" "$ROOT/scripts/prepare.py"
cp /usr/xtmp/lz280/packstar_flex_frontier38_20260920/prep_12651495/scripts/execute.py "$ROOT/scripts/execute.py"
cp "$REPO/slurm/scripts/run_markstar_backfill25_20260922.slurm" "$ROOT/scripts/run.slurm"
cp "${BASH_SOURCE[0]}" "$ROOT/scripts/prepare.slurm"
bash -n "$ROOT/scripts/run.slurm"
python3 "$ROOT/scripts/prepare.py" "$ROOT"
