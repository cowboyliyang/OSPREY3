#!/bin/bash
#SBATCH --job-name=flatsum_vs_markstar
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/flatsum_vs_markstar_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/flatsum_vs_markstar_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "FLAT_SUM BranchMARK* vs Original MARK*"
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

# ===== Run test =====
echo ""
echo "=========================================="
echo "Running testFlatSumVsOriginalMARKStar"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testFlatSumVsOriginalMARKStar" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx28g" 2>&1

# ===== Extract key results =====
echo ""
echo "=========================================="
echo "Key results (from XML report):"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestBranchMARKStar.xml
if [ -f "$XML" ]; then
    grep -oP 'time="[^"]*"' "$XML" | head -5
    grep -E 'TIMING|runtime:|Speedup|COMPARISON|MARK\*|BranchMARK|Branch decomposition|Branchwidth|TESS|falling back|epsilon|FLAT_SUM|Initialized|Finished' "$XML"
else
    echo "No report found"
fi

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
