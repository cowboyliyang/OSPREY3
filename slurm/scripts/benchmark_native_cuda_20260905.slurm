#!/bin/bash
#SBATCH --job-name=bench_native_cuda
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/bench_native_cuda_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/bench_native_cuda_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
set -euo pipefail

cd /home/users/lz280/OSPREY3-fresh-packstar
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$JAVA_HOME/bin:$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
export NVCC="$(readlink -f "$cuda_env/bin/nvcc")"
cuda_root="$(dirname "$(dirname "$NVCC")")"
run_root=/usr/xtmp/lz280/packstar_migration/bench-native-cuda-$SLURM_JOB_ID
mkdir -p "$run_root/tmp"
export TMPDIR="$run_root/tmp"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$run_root/gradle-build"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/native-gradle-home
export LD_LIBRARY_PATH=/usr/xtmp/lz280/packstar_migration/native-verify-12505273/lib:$cuda_root/lib64:$cuda_env/lib:${LD_LIBRARY_PATH:-}

test -f /usr/xtmp/lz280/packstar_migration/native-verify-12505273/lib/libCudaConfEcalc.so
./gradlew test --no-daemon --max-workers=1 \
    --project-cache-dir "$run_root/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    -Djna.library.path=/usr/xtmp/lz280/packstar_migration/native-verify-12505273/lib \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.benchmark_native_cuda_f64_throughput
