#!/bin/bash
#SBATCH --job-name=hybrid_iter
#SBATCH --output=hybrid_iter_%j.out
#SBATCH --error=hybrid_iter_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "Hybrid vs CCD Iteration Comparison"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean Gradle locks/daemons
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

# Compile
echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi

# Run continuous DP benchmark
echo ""
echo "=========================================="
echo "Running benchmarkHybridIterations"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmarkHybridIterations" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx24g" 2>&1

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
