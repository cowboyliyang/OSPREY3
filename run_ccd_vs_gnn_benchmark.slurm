#!/bin/bash
#SBATCH --job-name=ccd_vs_gnn
#SBATCH --output=ccd_vs_gnn_%j.out
#SBATCH --error=ccd_vs_gnn_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "CCD vs GNN Benchmark: 1CC8 8-res"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

# Isolate gradle state
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

NUM_CONFS="${NUM_CONFS:-100}"
TOP_K="${TOP_K:-10}"
MODEL_PATH="${MODEL_PATH:-gnn_data/1CC8_8res/model/gnn_model.onnx}"

echo "Config:"
echo "  numConfs=$NUM_CONFS"
echo "  topK=$TOP_K"
echo "  modelPath=$MODEL_PATH"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.eval.numConfs=$NUM_CONFS \
-Dosprey.gnn.eval.topK=$TOP_K \
-Dosprey.gnn.eval.modelPath=$MODEL_PATH"

echo ""
echo "Compiling..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmarkCCDvsGNN1CC8_8res" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx48g"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
