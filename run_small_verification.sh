#!/bin/bash
#SBATCH --job-name=dp_verify_small
#SBATCH --output=dp_verify_small_%j.log
#SBATCH --error=dp_verify_small_%j.err
#SBATCH --time=30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP Enhancements Verification (Small Scale)"
echo "2 flexible residues - Quick test"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 8GB, CPUs: 4, Time limit: 30min"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Stop any existing gradle daemons in this job
./gradlew --stop || true
sleep 3

echo "Compiling with DP enhancements..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful"
echo ""

# Run small-scale verification test
echo "Running DP enhancements verification test..."
echo "Configuration: 2 flexible residues (G649, A172)"
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"

TEST_RESULT=$?

echo ""
echo "End time: $(date)"

if [ $TEST_RESULT -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓✓✓ VERIFICATION SUCCESSFUL ✓✓✓"
    echo "=========================================="
    echo ""
    echo "Verified features:"
    echo "  ✓ Backtracking (getSelectedCorrectionsFromDP)"
    echo "  ✓ Statistics tracking (corrections considered/selected)"
    echo "  ✓ Enhanced printDPStats() output"
    echo ""
else
    echo ""
    echo "=========================================="
    echo "❌ VERIFICATION FAILED"
    echo "=========================================="
    exit 1
fi
