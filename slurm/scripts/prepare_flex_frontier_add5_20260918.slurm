#!/bin/bash
#SBATCH --job-name=ms_add5_prepare
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/packstar_flex_frontier12_20260918/prep_${SLURM_JOB_ID}
mkdir -p "$ROOT/tmp" "$ROOT/scripts"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
cp "$REPO/slurm/scripts/prepare_flex_frontier_add5_20260918.py" "$ROOT/scripts/prepare.py"
cp "$REPO/slurm/scripts/run_markstar_flex_add5_20260918.slurm" "$ROOT/scripts/run_markstar.slurm"
cp /usr/xtmp/lz280/packstar_flex_frontier10_20260917/prep_12631844/scripts/execute.py "$ROOT/scripts/execute.py"
cp "${BASH_SOURCE[0]}" "$ROOT/scripts/prepare.slurm"
python3 "$ROOT/scripts/prepare.py" "$ROOT"
