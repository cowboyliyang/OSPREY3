#!/bin/bash
#SBATCH --job-name=phase6_branchmarkstar
#SBATCH --output=phase6_branchmarkstar_%j.out
#SBATCH --error=phase6_branchmarkstar_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "Phase 6: BranchMARK* Correctness Test"
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
# Remove project-level lock files
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

# ===== Run both tests in a single Gradle invocation (preserves report) =====
echo ""
echo "=========================================="
echo "Running testBBKStarComparison_2RL0"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testBBKStarComparison_2RL0" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx28g" 2>&1

# ===== Extract key results from XML report (avoid dumping 250MB) =====
echo ""
echo "=========================================="
echo "Key results (from XML report):"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestBranchMARKStar.xml
if [ -f "$XML" ]; then
    grep -oP 'time="[^"]*"' "$XML" | head -5
    grep -E 'TIMING|runtime:|Speedup|CORRECTNESS|COMPARISON|FALLBACK|Standard-MARK|BranchMARK|sequences=|Matched:|Dense graph|Branch decomposition computed|Branchwidth|TESS|falling back|DP bounds|epsilon' "$XML"
else
    echo "No report found"
fi

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
