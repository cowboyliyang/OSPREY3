#!/bin/bash
#SBATCH --job-name=packstar-deadline-focus
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --gres=gpu:rtx_pro_6000:4
#SBATCH --exclusive
#SBATCH --mem=1000G
#SBATCH --time=06:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/packstar_deadline_focus_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/packstar_deadline_focus_%j.err
set -euo pipefail
: "${EXTENSION_PREP:?Set validated extension preparation path}"
export FOCUS_PACK_INDICES="${FOCUS_PACK_INDICES:-0,1,2,3,4,5,6,7,8,9,10,11}"
export RUN_ROOT=/usr/xtmp/lz280/markstar_deadline_focus_20261006/pack_${SLURM_JOB_ID}
mkdir -p "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export PACKSTAR_EXPECT_GPU='RTX PRO 6000'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/users/lz280/miniconda3/envs/confdiff/lib:${LD_LIBRARY_PATH:-}"
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
cp "${BASH_SOURCE[0]}" "$RUN_ROOT/submitted.slurm"
cp "$EXTENSION_PREP/run_group.py" "$RUN_ROOT/run_group.py"
lscpu > "$RUN_ROOT/lscpu.txt"
python3 "$RUN_ROOT/run_group.py"
