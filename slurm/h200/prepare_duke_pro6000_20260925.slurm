#!/bin/bash
#SBATCH --job-name=ps-pro6000-duke
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodelist=compsci-cluster-fitz-49
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:rtx_pro_6000:2
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_duke_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_duke_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588
ROOT=/usr/xtmp/lz280/packstar_pro6000_20260925/duke_${SLURM_JOB_ID}
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export TMPDIR="$ROOT/tmp" CUDA_CACHE_PATH="$ROOT/cuda_cache" PYTHONPYCACHEPREFIX="$ROOT/pycache"
mkdir -p "$ROOT/source/slurm/h200" "$ROOT/gate" "$ROOT/tmp" "$ROOT/resources/gpuKernels/cuda"
test -f "$BASE/READY"
# Reuse Duke's validated runtime. Only the GPU kernel resource is overlaid.
cp "$BASE/test_classpath.txt" "$ROOT/base_classpath.txt"
cp "$BASE/git_head.txt" "$ROOT/git_head.txt"
cp "$BASE/source.sha256" "$ROOT/base-source.sha256"
cp "$REPO/slurm/h200/"{run.py,frontier_active.tsv,baseline38.csv,Pro6000GpuGate.java} "$ROOT/source/slurm/h200/"
cp "$BASE/source/slurm/h200/production.properties" "$ROOT/source/slurm/h200/"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
cp "$REPO/slurm/h200/run_pro6000_20260925.slurm" "$ROOT/source/slurm/h200/"
bash -n "$ROOT/source/slurm/h200/run_pro6000_20260925.slurm"
CP="$(cat "$BASE/test_classpath.txt")"
"$JAVA_HOME/bin/javac" -cp "$CP" -d "$ROOT/gate" "$ROOT/source/slurm/h200/Pro6000GpuGate.java"
nvidia-smi > "$ROOT/gpus.txt"
# Test the existing Duke binary first; only recompile the exact same dp.cu if needed.
if "$JAVA_HOME/bin/java" --add-modules jdk.incubator.foreign -Xmx16g -Djava.io.tmpdir="$TMPDIR" \
    -cp "$ROOT/gate:$CP" Pro6000GpuGate > "$ROOT/existing_gpu_gate.log" 2>&1; then
    printf '%s\n' "$CP" > "$ROOT/test_classpath.txt"
else
    source /etc/profile.d/modules.sh 2>/dev/null || true
    if command -v module >/dev/null; then
        module load CUDA/12.8 2>/dev/null || module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
    fi
    unset NVCC CUDA_HOME
    for candidate in /usr/local/cuda/bin/nvcc /usr/local/cuda-12.8/bin/nvcc /usr/local/cuda-12.4/bin/nvcc; do
        if [ -x "$candidate" ]; then NVCC="$candidate"; break; fi
    done
    : "${NVCC:=$(command -v nvcc || true)}"
    test -x "$NVCC"
    "$NVCC" --version > "$ROOT/nvcc.txt"
    cp "$BASE/source/src/main/resources/gpuKernels/cuda/dp.cu" "$ROOT/resources/gpuKernels/cuda/dp.cu"
    # Same portable PTX target used by the repository's system-CUDA build path.
    "$NVCC" -fatbin -gencode=arch=compute_90,code=compute_90 \
        "$ROOT/resources/gpuKernels/cuda/dp.cu" -o "$ROOT/resources/gpuKernels/cuda/dp.bin"
    printf '%s:%s\n' "$ROOT/resources" "$CP" > "$ROOT/test_classpath.txt"
    "$JAVA_HOME/bin/java" --add-modules jdk.incubator.foreign -Xmx16g -Djava.io.tmpdir="$TMPDIR" \
        -cp "$ROOT/gate:$(cat "$ROOT/test_classpath.txt")" Pro6000GpuGate > "$ROOT/overlay_gpu_gate.log" 2>&1
fi
cd "$ROOT"
find source resources -type f -print0 | sort -z | xargs -0 sha256sum > source.sha256
printf 'base_build\t%s\nprepare_job\t%s\n' "$BASE" "$SLURM_JOB_ID" > build_manifest.tsv
touch READY
echo "READY=$ROOT"
