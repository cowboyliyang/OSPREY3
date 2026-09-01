#!/bin/bash
# Compile the migrated latest PackStar source tree under Slurm.
# All Gradle outputs, caches, and logs are placed under /usr/xtmp/lz280.
#SBATCH --job-name=packstar_migration_compile
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/slurm_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/gradle-build
GRADLE_ROOT=/usr/xtmp/lz280/packstar_migration/gradle-home

mkdir -p /usr/xtmp/lz280/packstar_migration
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT"
export GRADLE_USER_HOME="$GRADLE_ROOT"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load CUDA/12.8 >/dev/null 2>&1 \
        || module load CUDA/12.4 >/dev/null 2>&1 \
        || module load cuda/12.8 >/dev/null 2>&1 \
        || module load cuda/12.4 >/dev/null 2>&1 \
        || module load cuda/cuda-12.1 >/dev/null 2>&1 \
        || module load cuda >/dev/null 2>&1 \
        || true
fi
# The login environment may export a conda CUDA toolkit whose nvcc predates
# SM90.  Let the loaded system CUDA module (or the fixed system paths in the
# Gradle CUDA task) choose the compiler instead.
unset NVCC CUDA_HOME
echo "nvcc=$(command -v nvcc || echo missing)"
nvcc --version 2>&1 | tail -n 1 || true
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

cd "$REPO"
./gradlew compileCuda_dp compileJava \
    --no-daemon \
    --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g
