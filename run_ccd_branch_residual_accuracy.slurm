#!/bin/bash
#SBATCH --job-name=ccd_branch_acc
#SBATCH --output=ccd_branch_acc_%j.out
#SBATCH --error=ccd_branch_acc_%j.err
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
echo "Mode: branchResidual surrogate vs CCD/forcefield final eval"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

# isolate gradle state to avoid lock contention
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

# -------- Surrogate config (branchResidual) --------
TASK_TAG="benchmark1cc8_branch_residual_${SLURM_JOB_ID}"
CACHE_ROOT="/home/users/lz280/osprey-surrogate-cache"

# You can override these by exporting env vars before sbatch.
PRETRAIN_SAMPLES="${PRETRAIN_SAMPLES:-20000}"
FINETUNE_SAMPLES="${FINETUNE_SAMPLES:-2000}"
HASH_DIM="${HASH_DIM:-65536}"
LAMBDA_PRE="${LAMBDA_PRE:-0.01}"
LAMBDA_FINE="${LAMBDA_FINE:-0.1}"
CCD_NUM_CONFS="${CCD_NUM_CONFS:-15}"
CCD_MAX_ITER="${CCD_MAX_ITER:-30}"
CCD_TOPK="${CCD_TOPK:-5}"

echo "Surrogate config:"
echo "  enabled=true"
echo "  type=branchResidual"
echo "  cacheRoot=$CACHE_ROOT"
echo "  taskTag=$TASK_TAG"
echo "  branchResidual.enabled=true"
echo "  branchResidual.samplerMode=sparse"
echo "  branchResidual.pretrain.griddpSamples=$PRETRAIN_SAMPLES"
echo "  branchResidual.finetune.ccdSamples=$FINETUNE_SAMPLES"
echo "  branchResidual.ridge.lambdaPre=$LAMBDA_PRE"
echo "  branchResidual.ridge.lambdaFine=$LAMBDA_FINE"
echo "  branchResidual.feature.hashDim=$HASH_DIM"
echo "  ccdOnly.numConfs=$CCD_NUM_CONFS"
echo "  ccdOnly.maxIter=$CCD_MAX_ITER"
echo "  ccdOnly.topK=$CCD_TOPK"

export JAVA_TOOL_OPTIONS="\
-Dosprey.bench.surrogate.enabled=true \
-Dosprey.bench.surrogate.type=branchResidual \
-Dosprey.bench.surrogate.cacheRoot=$CACHE_ROOT \
-Dosprey.bench.surrogate.taskTag=$TASK_TAG \
-Dosprey.bench.surrogate.forceRetrain=false \
-Dosprey.bench.surrogate.branchResidual.enabled=true \
-Dosprey.bench.surrogate.branchResidual.samplerMode=sparse \
-Dosprey.bench.surrogate.branchResidual.pretrain.griddpSamples=$PRETRAIN_SAMPLES \
-Dosprey.bench.surrogate.branchResidual.finetune.ccdSamples=$FINETUNE_SAMPLES \
-Dosprey.bench.surrogate.branchResidual.ridge.lambdaPre=$LAMBDA_PRE \
-Dosprey.bench.surrogate.branchResidual.ridge.lambdaFine=$LAMBDA_FINE \
-Dosprey.bench.surrogate.branchResidual.feature.hashDim=$HASH_DIM \
-Dosprey.bench.surrogate.branchResidual.pretrain.useFinalEval=true \
-Dosprey.bench.surrogate.branchResidual.gridSize=2 \
-Dosprey.bench.surrogate.branchResidual.gridThreads=1 \
-Dosprey.bench.surrogate.branchResidual.lowEnergyPoolLimit=4096 \
-Dosprey.bench.surrogate.branchResidual.randomSeed=12345 \
-Dosprey.bench.surrogate.branchResidual.ridge.maxIter=256 \
-Dosprey.bench.surrogate.branchResidual.ridge.tol=1e-6 \
-Dosprey.bench.surrogate.ccdOnly.numConfs=$CCD_NUM_CONFS \
-Dosprey.bench.surrogate.ccdOnly.maxIter=$CCD_MAX_ITER \
-Dosprey.bench.surrogate.ccdOnly.topK=$CCD_TOPK"

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
