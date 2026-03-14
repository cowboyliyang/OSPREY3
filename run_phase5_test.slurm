#!/bin/bash
#SBATCH --job-name=phase5_test
#SBATCH --output=phase5_test_%j.out
#SBATCH --error=phase5_test_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "Phase 5: Partial Start Cache Performance Test"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean Gradle daemons and caches
echo "Cleaning Gradle..."
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true

# Compile
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon

# Run test (--rerun-tasks forces re-execution even if cached)
echo ""
echo "Running Phase 5 test..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestPartialStartCache.testPartialStartCacheIsolated" --no-daemon --rerun-tasks

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
