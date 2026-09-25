#!/bin/bash
#SBATCH --job-name=verify_native_ccd
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/verify_native_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/verify_native_%j.err
#SBATCH --time=02:00:00
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
# Invoke nvcc inside its toolkit: the conda bin symlink otherwise makes nvcc
# search conda/lib64, while the complete toolkit libraries live below pkgs.
export NVCC="$(readlink -f "$cuda_env/bin/nvcc")"
cuda_root="$(dirname "$(dirname "$NVCC")")"
export PATH="$cuda_root/bin:$PATH"
run_root=/usr/xtmp/lz280/packstar_migration/native-verify-$SLURM_JOB_ID
mkdir -p "$run_root/tmp" "$run_root/lib"
export TMPDIR="$run_root/tmp"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="${BUILD_ROOT:-$run_root/gradle-build}"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/native-gradle-home
native_lib_dir=${NATIVE_LIBRARY_DIR:-$run_root/lib}
export LD_LIBRARY_PATH="$native_lib_dir:$cuda_root/lib64:$cuda_env/lib:${LD_LIBRARY_PATH:-}"
"$NVCC" --version
if [[ -z "${NATIVE_LIBRARY_DIR:-}" ]]; then
cmake -S src/main/cu/CudaConfEcalc -B "$run_root/cmake" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER="$NVCC" \
    -DCMAKE_CUDA_ARCHITECTURES=OFF -DLIBRARY_OUTPUT_PATH="$run_root/lib"
cmake --build "$run_root/cmake" --target CudaConfEcalc -j1
fi
test -f "$native_lib_dir/libCudaConfEcalc.so"
./gradlew test --no-daemon --max-workers=1 \
    --project-cache-dir "${PROJECT_CACHE:-$run_root/project-cache}" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    -Djna.library.path="$native_lib_dir" \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimize_cuda_2RL0_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimize_cuda_1DG9_6f_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimizeEnergies_cuda_2RL0_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimizeEnergies_cuda_1DG9_6f_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.audit_cuda_cpu_f64_consistency
