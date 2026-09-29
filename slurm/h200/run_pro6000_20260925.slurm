#!/bin/bash
#SBATCH --job-name=ps-pro6000-active
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:rtx_pro_6000:4
#SBATCH --time=01:35:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_active_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_active_%A_%a.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${SLURM_ARRAY_TASK_ID:?Supply design index}"
: "${BUILD_ROOT:?Validated PRO 6000 build}"
export REPO=/home/users/lz280/OSPREY3-fresh-packstar
export INPUT_ROOT=/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package
export RESULT_ROOT=/usr/xtmp/lz280/packstar_pro6000_20260925/results
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PACKSTAR_EXPECT_GPU='RTX PRO 6000'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/users/lz280/miniconda3/envs/confdiff/lib:${LD_LIBRARY_PATH:-}"
test -f "$BUILD_ROOT/READY"
DRIVER_ROOT="$RESULT_ROOT/drivers/J${SLURM_JOB_ID}_T${SLURM_ARRAY_TASK_ID}"
mkdir -p "$RESULT_ROOT/drivers"
mkdir "$DRIVER_ROOT"
cp "$BUILD_ROOT/source/slurm/h200/run.py" "$DRIVER_ROOT/run.py"
cp "$BUILD_ROOT/source/slurm/h200/frontier_active.tsv" "$DRIVER_ROOT/frontier_active.tsv"
cp "${BASH_SOURCE[0]}" "$DRIVER_ROOT/submitted.slurm"
export PYTHONPYCACHEPREFIX="$DRIVER_ROOT/pycache"
export TMPDIR="$DRIVER_ROOT/tmp"
mkdir "$TMPDIR"
cd "$INPUT_ROOT"
sha256sum -c SHA256SUMS > "$DRIVER_ROOT/input_checksum.txt"
# Every design must match the frozen positions and sequence order before searching.
timeout 30m python3 "$DRIVER_ROOT/run.py" --frontier-manifest "$DRIVER_ROOT/frontier_active.tsv" \
    --mode preflight --arm pair-only --rb 1 --heap-gib 8 --host-gib 4 --gpus 4 --gpu-gib 85
# The node worker enforces a combined 90-minute limit for this entire script.
python3 "$DRIVER_ROOT/run.py" \
    --frontier-manifest "$DRIVER_ROOT/frontier_active.tsv" \
    --arm pair-only --rb 1 --gpus 4 --gpu-gib 85 --heap-gib 850 --host-gib 800
