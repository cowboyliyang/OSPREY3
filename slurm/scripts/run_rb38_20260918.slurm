#!/bin/bash
#SBATCH --job-name=packstar_rb38
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=112G
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -f "${SWEEP_ROOT:?}/READY"
export BUILD_ROOT=/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056
export INPUT_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export RESULT_ROOT="$SWEEP_ROOT/rb${RB:?}"
export PACKSTAR_EXPECT_GPU=''
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
# Budget 80% of the one allocated GPU's memory; never inspect other jobs' GPUs.
gpu_mib=$(nvidia-smi -i "${CUDA_VISIBLE_DEVICES:?}" --query-gpu=memory.total --format=csv,noheader,nounits)
gpu_gib=$((gpu_mib * 80 / 100 / 1024))
test "$gpu_gib" -ge 1
exec python3 "$SWEEP_ROOT/source/run.py" --cohort baseline38 --mode full \
    --arm budget-forward --seed 42 --rb "$RB" --gpus 1 \
    --heap-gib 96 --host-gib 80 --gpu-gib "$gpu_gib"
