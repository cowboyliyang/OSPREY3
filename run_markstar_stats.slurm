#!/bin/bash
#SBATCH --job-name=markstar_stats
#SBATCH --output=markstar_stats_%j.out
#SBATCH --error=markstar_stats_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "MARK* Leaf vs Internal Statistics Test"
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

# ===== Run MARK* test =====
echo ""
echo "=========================================="
echo "Running testMARKStarVsKStar (MARK* only, 8 flex, epsilon=0.68)"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestMARKStar.testMARKStarVsKStar" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx12g" 2>&1

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
