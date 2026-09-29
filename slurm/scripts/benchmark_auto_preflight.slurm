#!/bin/bash
# Paired full-design run on an exclusive eight-A5000 node. Existing EMAT input ~few MiB/design.
#SBATCH --job-name=ps-auto-benefit
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=fennario-02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --gres=gpu:a5000:8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=04:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
export AUTO_BENCH_ROOT=/usr/xtmp/lz280/packstar_auto_preflight/benefit_${SLURM_JOB_ID}
export AUTO_PREFLIGHT_CLASSPATH=${AUTO_PREFLIGHT_CLASSPATH:-/usr/xtmp/lz280/packstar_auto_preflight/query_fix_12705178/test_classpath.txt}
export AUTO_BENCH_SOURCE=/usr/xtmp/lz280/packstar_auto_preflight/check_12705173/source
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="/home/users/lz280/miniconda3/envs/confdiff/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR="$AUTO_BENCH_ROOT/tmp"
export CUDA_CACHE_PATH="$AUTO_BENCH_ROOT/cuda_cache"
export PYTHONPYCACHEPREFIX="$AUTO_BENCH_ROOT/pycache"
mkdir -p "$TMPDIR" "$CUDA_CACHE_PATH" "$AUTO_BENCH_ROOT/gate"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/benchmark_auto_preflight.py "$AUTO_BENCH_ROOT/driver.py"
cp "${BASH_SOURCE[0]}" "$AUTO_BENCH_ROOT/submitted.slurm"
nvidia-smi -i "$CUDA_VISIBLE_DEVICES" > "$AUTO_BENCH_ROOT/gpus.txt"
lscpu > "$AUTO_BENCH_ROOT/cpus.txt"
CP=$(<"$AUTO_PREFLIGHT_CLASSPATH")
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/Pro6000GpuGate.java "$AUTO_BENCH_ROOT/gate/"
javac -cp "$CP" -d "$AUTO_BENCH_ROOT/gate" "$AUTO_BENCH_ROOT/gate/Pro6000GpuGate.java"
java --add-modules jdk.incubator.foreign -Djava.io.tmpdir="$TMPDIR" \
    -cp "$AUTO_BENCH_ROOT/gate:$CP" Pro6000GpuGate > "$AUTO_BENCH_ROOT/gpu_gate.log" 2>&1
python3 "$AUTO_BENCH_ROOT/driver.py"
