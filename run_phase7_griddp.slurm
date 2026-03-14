#!/bin/bash
#SBATCH --job-name=phase7_griddp
#SBATCH --output=phase7_griddp_%j.out
#SBATCH --error=phase7_griddp_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "Phase 7: Grid DP + MARK* Test"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean Gradle
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

# Run Phase 7 test
echo ""
echo "=========================================="
echo "Running testGridDPPhase7"
echo "=========================================="
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testGridDPPhase7" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx12g" 2>&1

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
