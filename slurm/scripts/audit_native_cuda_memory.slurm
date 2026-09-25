#!/bin/bash
# Slurm-only native CUDA regression/audit.  This exercises both single-job and
# native batched minimization after the memory-aware stream sizing changes.
#SBATCH --job-name=packstar_native_cuda_audit
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/native_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/native_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/native-gradle-build
GRADLE_ROOT=/usr/xtmp/lz280/packstar_migration/native-gradle-home
NATIVE_BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/native-cuda-build
NATIVE_LIB=/usr/xtmp/lz280/packstar_migration/libCudaConfEcalc.so
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
CUDA_ENV=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$JAVA_HOME/bin:$CUDA_ENV/bin:$CUDA_ENV/cuda/bin:$PATH"
export NVCC="$CUDA_ENV/bin/nvcc"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT"
export GRADLE_USER_HOME="$GRADLE_ROOT"
export LD_LIBRARY_PATH="$REPO/src/main/resources/linux-x86-64:$CUDA_ENV/lib:$CUDA_ENV/cuda/lib64:${LD_LIBRARY_PATH:-}"

cd "$REPO"

# Build the native CUDA library in the Slurm scratch area.  CMake's historical
# project file writes the shared object into the resource directory; move that
# one explicitly into xtmp and leave a stable resource-side symlink for JNA.
cmake -S src/main/cu/CudaConfEcalc -B "$NATIVE_BUILD_ROOT" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_COMPILER="$NVCC" \
    -DCMAKE_CUDA_ARCHITECTURES=OFF
cmake --build "$NATIVE_BUILD_ROOT" --target CudaConfEcalc -j1

source_lib="$REPO/src/main/resources/linux-x86-64/libCudaConfEcalc.so"
test -f "$source_lib"
if test -L "$source_lib"; then
    resolved=$(readlink -f "$source_lib")
    test "$resolved" = "$NATIVE_LIB"
else
    cp -f "$source_lib" "$NATIVE_LIB"
    cmp -s "$source_lib" "$NATIVE_LIB"
    mv "$source_lib" "$NATIVE_LIB"
    ln -s "$NATIVE_LIB" "$source_lib"
fi

./gradlew test --no-daemon --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimize_cuda_2RL0_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimize_cuda_1DG9_6f_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimizeEnergies_cuda_2RL0_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.minimizeEnergies_cuda_1DG9_6f_f64 \
    --tests edu.duke.cs.osprey.energy.compiled.TestNativeConfEnergyCalculator.audit_cuda_cpu_f64_consistency
