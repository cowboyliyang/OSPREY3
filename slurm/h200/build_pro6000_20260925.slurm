#!/bin/bash
#SBATCH --job-name=ps-pro6000-build
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodelist=compsci-cluster-fitz-48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:rtx_pro_6000:1
#SBATCH --time=02:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_build_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_build_%j.err
set -euo pipefail
export REPO=/home/users/lz280/OSPREY3-fresh-packstar
export BUILD_ROOT=/usr/xtmp/lz280/packstar_pro6000_20260925/build_${SLURM_JOB_ID}
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_CUDA_ARCHS=120
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export NVCC="$(readlink -f "$cuda_env/bin/nvcc")"
export PATH="$JAVA_HOME/bin:$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
export LD_LIBRARY_PATH="$cuda_env/lib:${LD_LIBRARY_PATH:-}"
test -x "$NVCC"
bash "$REPO/slurm/h200/build.slurm"
