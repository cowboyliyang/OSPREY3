#!/bin/bash
#SBATCH --job-name=bench_pro6000
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/bench_pro6000_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/bench_pro6000_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# The cluster exposes the Pro 6000 nodes in compsci-gpu under this GRES name.
#SBATCH --partition=compsci-gpu
#SBATCH --account=grisman
#SBATCH --gres=gpu:rtx_pro_6000:1
set -euo pipefail

cd /home/users/lz280/OSPREY3-fresh-packstar
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$JAVA_HOME/bin:$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
export NVCC="$(readlink -f "$cuda_env/bin/nvcc")"
cuda_root="$(dirname "$(dirname "$NVCC")")"
run_root=/usr/xtmp/lz280/packstar_migration/bench-pro6000-$SLURM_JOB_ID
mkdir -p "$run_root/tmp" "$run_root/lib"
export TMPDIR="$run_root/tmp"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$run_root/gradle-build"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/native-gradle-home
export LD_LIBRARY_PATH="$run_root/lib:$cuda_root/lib64:$cuda_env/lib:${LD_LIBRARY_PATH:-}"

nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader
"$NVCC" --version
cmake -S src/main/cu/CudaConfEcalc -B "$run_root/cmake" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER="$NVCC" \
    -DCMAKE_CUDA_ARCHITECTURES=OFF -DLIBRARY_OUTPUT_PATH="$run_root/lib"
cmake --build "$run_root/cmake" --target CudaConfEcalc -j1
test -f "$run_root/lib/libCudaConfEcalc.so"

./gradlew test --no-daemon --max-workers=1 \
    --project-cache-dir "$run_root/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    -Djna.library.path="$run_root/lib" \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.benchmark_native_cuda_f64_throughput
