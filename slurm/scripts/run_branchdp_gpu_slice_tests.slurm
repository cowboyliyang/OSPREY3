#!/bin/bash
#SBATCH --job-name=branchdp_gpu_full_matrix
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_full_matrix_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_full_matrix_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

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
    module load cuda/12.8 >/dev/null 2>&1 || module load cuda/12.4 >/dev/null 2>&1 || module load cuda >/dev/null 2>&1 || true
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
    if command -v nvcc >/dev/null 2>&1; then candidates+=("$(command -v nvcc)"); fi

    local nvcc home
    for nvcc in "${candidates[@]}"; do
        if [ ! -x "$nvcc" ]; then continue; fi
        home="$(cd "$(dirname "$nvcc")/.." && pwd)"
        if [ -f "$home/include/cuda_runtime.h" ]; then
            export NVCC="$nvcc"
            export CUDA_HOME="$home"
            return 0
        fi
    done
    return 1
}

if ! pick_nvcc; then
    echo "Could not find a complete CUDA toolkit with cuda_runtime.h" >&2
    echo "NVCC=${NVCC:-unset}" >&2
    echo "CUDA_HOME=${CUDA_HOME:-unset}" >&2
    command -v nvcc || true
    exit 1
fi

cuda_path_entries=()
for d in "$CUDA_HOME/bin" "$CUDA_HOME/nvvm/bin" "$CUDA_HOME/../pkgs/cuda-toolkit/bin" "$CUDA_HOME/../pkgs/cuda-toolkit/nvvm/bin"; do
    if [ -d "$d" ]; then
        cuda_path_entries+=("$(cd "$d" && pwd)")
    fi
done
if [ "${#cuda_path_entries[@]}" -gt 0 ]; then
    export PATH="$(IFS=:; echo "${cuda_path_entries[*]}"):$PATH"
fi

cuda_include_entries=()
for d in "$CUDA_HOME/include" "$CUDA_HOME/targets/x86_64-linux/include"; do
    if [ -d "$d" ]; then
        cuda_include_entries+=("$(cd "$d" && pwd)")
    fi
done
if [ "${#cuda_include_entries[@]}" -gt 0 ]; then
    export CPATH="$(IFS=:; echo "${cuda_include_entries[*]}")${CPATH:+:$CPATH}"
fi

echo "=============================================="
echo "Branch-DP two-GPU full regression matrix"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "Start: $(date)"
echo "Java: $(java -version 2>&1 | head -n 1)"
echo "NVCC: $NVCC"
"$NVCC" --version | tail -n 1
echo "CUDA_HOME: $CUDA_HOME"
echo "CUDA PATH entries: ${cuda_path_entries[*]:-none}"
echo "CUDA include entries: ${cuda_include_entries[*]:-none}"
echo "GRADLE_USER_HOME: $GRADLE_USER_HOME"
echo "=============================================="

nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader

./gradlew compileCuda_dp \
    --offline \
    --no-daemon \
    --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g

./gradlew compileTestJava \
    --offline \
    --no-daemon \
    --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g

./gradlew test \
    --offline \
    --no-daemon \
    --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g \
    -DtestMaxHeap=16g \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP \
    --tests edu.duke.cs.osprey.branchdp.TestComputeFullDPSharded \
    --tests edu.duke.cs.osprey.branchdp.TestDPTable

echo ""
echo "End: $(date)"
