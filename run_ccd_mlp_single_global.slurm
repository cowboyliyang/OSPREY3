#!/bin/bash
#SBATCH --job-name=ccd_mlp_single
#SBATCH --output=ccd_mlp_single_%j.out
#SBATCH --error=ccd_mlp_single_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

echo "=============================================="
echo "BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy"
echo "Mode: single-task-global MLP surrogate (CCD-only)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

# isolate gradle state to avoid lock contention with other jobs
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

TASK_TAG="benchmark1cc8_single_global_${SLURM_JOB_ID}"
CACHE_ROOT="/home/users/lz280/osprey-surrogate-cache"

echo "Surrogate config:"
echo "  enabled=true"
echo "  type=mlp"
echo "  mlp.singleTaskModel=true"
echo "  cacheRoot=$CACHE_ROOT"
echo "  taskTag=$TASK_TAG"
echo "  mlp.globalNumTuples=800"
echo "  mlp.globalSamplesPerTuple=48"
echo "  mlp.hidden1=128"
echo "  mlp.hidden2=128"
echo "  mlp.epochs=1000"
echo "  mlp.batchSize=256"
echo "  mlp.learningRate=0.001"
echo "  ccdOnly.numConfs=15"
echo "  ccdOnly.maxIter=30"
echo "  requireFullApproximation=true"
echo "  gridDPFallbackToForcefield=false"

export JAVA_TOOL_OPTIONS="\
-Dosprey.bench.surrogate.enabled=true \
-Dosprey.bench.surrogate.type=mlp \
-Dosprey.bench.surrogate.mlp.singleTaskModel=true \
-Dosprey.bench.surrogate.cacheRoot=$CACHE_ROOT \
-Dosprey.bench.surrogate.taskTag=$TASK_TAG \
-Dosprey.bench.surrogate.mlp.globalNumTuples=800 \
-Dosprey.bench.surrogate.mlp.globalSamplesPerTuple=48 \
-Dosprey.bench.surrogate.mlp.hidden1=128 \
-Dosprey.bench.surrogate.mlp.hidden2=128 \
-Dosprey.bench.surrogate.mlp.epochs=1000 \
-Dosprey.bench.surrogate.mlp.batchSize=256 \
-Dosprey.bench.surrogate.mlp.learningRate=0.001 \
-Dosprey.bench.surrogate.ccdOnly.numConfs=15 \
-Dosprey.bench.surrogate.ccdOnly.maxIter=30 \
-Dosprey.bench.surrogate.requireFullApproximation=true \
-Dosprey.bench.surrogate.gridDPFallbackToForcefield=false"

echo ""
echo "Compiling..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava
echo "Compilation successful."

echo ""
echo "=========================================="
echo "Running BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy"
echo "=========================================="
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx96g"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
