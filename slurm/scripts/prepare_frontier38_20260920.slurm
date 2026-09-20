#!/bin/bash
#SBATCH --job-name=ms38_prepare
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --exclude=fennario-[01-06]
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/packstar_flex_frontier38_20260920/prep_${SLURM_JOB_ID}
mkdir -p "$ROOT/tmp" "$ROOT/scripts"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
cp "$REPO/slurm/scripts/prepare_frontier38_20260920.py" "$ROOT/scripts/prepare.py"
cp "$REPO/slurm/scripts/execute_markstar_flex_add4_20260918.py" "$ROOT/scripts/execute.py"
cp "$REPO/slurm/scripts/run_markstar_frontier38_20260920.slurm" "$ROOT/scripts/run.slurm"
cp "${BASH_SOURCE[0]}" "$ROOT/scripts/prepare.slurm"
bash -n "$ROOT/scripts/run.slurm"
python3 "$ROOT/scripts/prepare.py" "$ROOT"
