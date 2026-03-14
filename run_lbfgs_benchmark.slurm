#!/bin/bash
#SBATCH --job-name=lbfgs_bench
#SBATCH --output=lbfgs_benchmark_%j.out
#SBATCH --error=lbfgs_benchmark_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "L-BFGS Hybrid Benchmark"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

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

# ===== Run benchmarkLBFGSHybrid =====
echo ""
echo "=========================================="
echo "Running benchmarkLBFGSHybrid"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmarkLBFGSHybrid" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx48g" 2>&1

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
