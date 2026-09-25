#!/bin/bash
#SBATCH --job-name=packstar_ccd_scalar_debug
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/scalar_ccd_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/scalar_ccd_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

set -euo pipefail
cd /home/users/lz280/OSPREY3-fresh-packstar

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
if [[ -n "${CUDA_MODULE:-}" ]]; then
    source /etc/profile.d/modules.sh
    module load "$CUDA_MODULE"
fi
export NVCC="${NVCC_OVERRIDE:-$cuda_env/bin/nvcc}"
export CUDA_LAUNCH_BLOCKING=1
export OSPREY_GRADLE_XTMP_BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/scalar-gradle-build
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/scalar-gradle-home

./gradlew compileCuda_residueCcd compileJava compileTestJava \
  --no-daemon --max-workers=1 \
  --init-script slurm/scripts/gradle_xtmp.init.gradle \
  -Dorg.gradle.jvmargs=-Xmx8g

./gradlew test \
  --tests edu.duke.cs.osprey.minimization.TestMinimization.testResidueCudaCCD1Stream \
  --tests edu.duke.cs.osprey.minimization.TestMinimization.testResidueCudaCCD2Stream \
  --no-daemon --max-workers=1 \
  --init-script slurm/scripts/gradle_xtmp.init.gradle \
  -Dorg.gradle.jvmargs=-Xmx8g --stacktrace
