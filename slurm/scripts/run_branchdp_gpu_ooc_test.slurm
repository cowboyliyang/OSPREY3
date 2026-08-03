#!/bin/bash
#SBATCH --job-name=branchdp_gpu_ooc
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_ooc_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_ooc_%j.err
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:2

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="${SLURM_TMPDIR:-/tmp}/osprey3-gradle-${SLURM_JOB_ID}"

mkdir -p "$GRADLE_USER_HOME/wrapper" "$GRADLE_USER_HOME/caches"
if [ -d "$HOME/.gradle/wrapper" ]; then
    cp -a "$HOME/.gradle/wrapper/." "$GRADLE_USER_HOME/wrapper/" 2>/dev/null || true
fi
if [ -d "$HOME/.gradle/caches/modules-2" ]; then
    cp -a "$HOME/.gradle/caches/modules-2" "$GRADLE_USER_HOME/caches/" 2>/dev/null || true
fi

if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load cuda/12.8 >/dev/null 2>&1 \
        || module load cuda/12.4 >/dev/null 2>&1 \
        || module load cuda/cuda-12.1 >/dev/null 2>&1 \
        || module load cuda >/dev/null 2>&1 \
        || true
fi

pick_nvcc() {
    local candidates=()
    if [ -n "${NVCC:-}" ]; then candidates+=("$NVCC"); fi
    if [ -n "${CUDA_HOME:-}" ]; then candidates+=("$CUDA_HOME/bin/nvcc"); fi
    candidates+=(
        /usr/local/cuda/bin/nvcc
        /usr/local/cuda-12.8/bin/nvcc
        /usr/local/cuda-12.4/bin/nvcc
        /usr/local/cuda-12.2/bin/nvcc
    )
    if command -v nvcc >/dev/null 2>&1; then
        candidates+=("$(command -v nvcc)")
    fi

    local candidate home
    for candidate in "${candidates[@]}"; do
        if [ ! -x "$candidate" ]; then continue; fi
        home="$(cd "$(dirname "$candidate")/.." && pwd)"
        if [ -f "$home/include/cuda_runtime.h" ] \
            || [ -f "$home/targets/x86_64-linux/include/cuda_runtime.h" ]; then
            export NVCC="$candidate"
            export CUDA_HOME="$home"
            return 0
        fi
    done
    return 1
}

if ! pick_nvcc; then
    echo "No nvcc found on $SLURM_NODELIST" >&2
    exit 1
fi

cuda_path_entries=()
for directory in "$CUDA_HOME/bin" "$CUDA_HOME/nvvm/bin" \
    "$CUDA_HOME/../pkgs/cuda-toolkit/bin" \
    "$CUDA_HOME/../pkgs/cuda-toolkit/nvvm/bin"; do
    if [ -d "$directory" ]; then
        cuda_path_entries+=("$(cd "$directory" && pwd)")
    fi
done
if [ "${#cuda_path_entries[@]}" -gt 0 ]; then
    export PATH="$(IFS=:; echo "${cuda_path_entries[*]}"):$PATH"
fi

cuda_includes=()
for directory in "$CUDA_HOME/include" "$CUDA_HOME/targets/x86_64-linux/include"; do
    if [ -d "$directory" ]; then
        cuda_includes+=("$(cd "$directory" && pwd)")
    fi
done
if [ "${#cuda_includes[@]}" -gt 0 ]; then
    export CPATH="$(IFS=:; echo "${cuda_includes[*]}")${CPATH:+:$CPATH}"
fi

echo "node=$SLURM_NODELIST gpus=${CUDA_VISIBLE_DEVICES:-unset} start=$(date)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
"$NVCC" --version | tail -n 1

./gradlew --project-cache-dir "$GRADLE_USER_HOME/project-cache" \
    compileCuda_dp compileTestJava --offline --no-daemon --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g

echo "git_head=$(git rev-parse HEAD)"
sha256sum \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuFullDP.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.java \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuFullDP.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.class \
    build/resources/main/gpuKernels/cuda/dp.bin

./gradlew --project-cache-dir "$GRADLE_USER_HOME/project-cache" \
    test --offline --no-daemon --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g -DtestMaxHeap=4g \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.outOfCorePlannerAndPackingBoundBothRowAndLambdaAxes_withoutGpu \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.outOfCoreStablePartialMergeHandlesExtremeNanAndNegativeInfinity_withoutGpu \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCoreExtremeNanNegativeInfinityAndAllocationAudit \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_bilateralRowAndLambdaTilingMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_fileBackedChildrenAndOutputMatchJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_multiGpuMatchesJava

echo "end=$(date)"
