#!/bin/bash
# Submit with explicit --partition, --gres, --nodelist, --job-name, and --output.
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
PROJECT_CACHE=${SLURM_TMPDIR:-/tmp}/osprey3-calibration-cache-${SLURM_JOB_ID}

cd "$REPO"
export JAVA_HOME
export PATH="$JAVA_HOME/bin:$PATH"

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

echo "calibration_job=$SLURM_JOB_ID node=$SLURMD_NODENAME gpus=${CUDA_VISIBLE_DEVICES:-unset} start=$(date --iso-8601=seconds)"
nvidia-smi --query-gpu=name,uuid,memory.total,memory.free \
    --format=csv,noheader
"$NVCC" --version | tail -n 1

./gradlew --project-cache-dir "$PROJECT_CACHE" \
    compileCuda_dp compileTestJava --offline --no-daemon --max-workers=4 \
    -Dorg.gradle.jvmargs=-Xmx3g

echo "git_head=$(git rev-parse HEAD)"
sha256sum \
    src/test/java/edu/duke/cs/osprey/branchdp/TestGpuFullDP.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuFullDP.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.java \
    build/classes/java/test/edu/duke/cs/osprey/branchdp/TestGpuFullDP.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuFullDP.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.class \
    build/resources/main/gpuKernels/cuda/dp.bin

./gradlew --project-cache-dir "$PROJECT_CACHE" \
    test --offline --no-daemon --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx3g -DtestMaxHeap=8g \
    -Dbranchdp.test.gpu.calibrate=true \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuLeafThroughputCalibration \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCoreThroughputCalibration

echo "calibration_job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
