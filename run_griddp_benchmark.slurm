#!/bin/bash
#SBATCH --job-name=ccd_mlp_accuracy
#SBATCH --output=ccd_mlp_accuracy_%j.out
#SBATCH --error=ccd_mlp_accuracy_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy"
echo "7-residue confspace, CCD-only MLP accuracy vs forcefield final eval"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH=$JAVA_HOME/bin:$PATH

TASK_TAG="benchmark1cc8_mlp_${SLURM_JOB_ID}"
CACHE_ROOT="/home/users/lz280/osprey-surrogate-cache"

echo "Surrogate config:"
echo "  enabled=true"
echo "  type=mlp"
echo "  cacheRoot=$CACHE_ROOT"
echo "  taskTag=$TASK_TAG"
echo "  samplesPerParam=16"
echo "  mlp.hidden1=64"
echo "  mlp.hidden2=64"
echo "  mlp.minSamplesPerModel=512"
echo "  mlp.maxSamplesPerModel=20000"
echo "  mlp.epochs=800"
echo "  mlp.batchSize=64"
echo "  mlp.learningRate=0.001"
echo "  requireFullApproximation=true"
echo "  gridDPFallbackToForcefield=false"

export JAVA_TOOL_OPTIONS="\
-Dosprey.bench.surrogate.enabled=true \
-Dosprey.bench.surrogate.type=mlp \
-Dosprey.bench.surrogate.cacheRoot=$CACHE_ROOT \
-Dosprey.bench.surrogate.taskTag=$TASK_TAG \
-Dosprey.bench.surrogate.samplesPerParam=16 \
-Dosprey.bench.surrogate.mlp.hidden1=64 \
-Dosprey.bench.surrogate.mlp.hidden2=64 \
-Dosprey.bench.surrogate.mlp.minSamplesPerModel=512 \
-Dosprey.bench.surrogate.mlp.maxSamplesPerModel=20000 \
-Dosprey.bench.surrogate.mlp.epochs=800 \
-Dosprey.bench.surrogate.mlp.batchSize=64 \
-Dosprey.bench.surrogate.mlp.learningRate=0.001 \
-Dosprey.bench.surrogate.ccdOnly.numConfs=15 \
-Dosprey.bench.surrogate.ccdOnly.maxIter=30 \
-Dosprey.bench.surrogate.requireFullApproximation=true \
-Dosprey.bench.surrogate.gridDPFallbackToForcefield=false"

# ===== Clean Gradle locks and daemons =====
echo ""
echo "Cleaning Gradle locks and daemons..."
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

# ===== Compile =====
echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation successful."

# ===== Run CCD-only MLP accuracy benchmark =====
echo ""
echo "=========================================="
echo "Running BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmark1CC8CCDOnlyMLPAccuracy" \
    --no-daemon \
    --max-workers=1 \
    -Dorg.gradle.jvmargs="-Xmx96g" 2>&1

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
