#!/bin/bash
#SBATCH --job-name=dp_phase1
#SBATCH --output=dp_phase1_%j.log
#SBATCH --error=dp_phase1_%j.err
#SBATCH --time=45:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP-MARKStar Phase 1 Optimizations Test"
echo "Binary Search + Parallel DP"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 16GB, CPUs: 4, Time limit: 45min"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Stop any existing gradle daemons
./gradlew --stop || true
sleep 3

echo "Compiling with Phase 1 optimizations..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful"
echo ""

echo "=========================================="
echo "Test: Compare Original DP vs Optimized DP"
echo "  - Original DP: O(n²) complexity"
echo "  - Optimized DP: O(n log n) with binary search"
echo "  - Parallel DP: Multi-threaded for large problems"
echo "=========================================="
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPOptimizations.testPhase1Optimizations"

TEST_RESULT=$?

echo ""
echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo "End time: $(date)"
echo ""

if [ $TEST_RESULT -eq 0 ]; then
    echo "✓✓✓ PHASE 1 OPTIMIZATIONS TEST PASSED ✓✓✓"
    echo ""
    echo "Optimizations Verified:"
    echo "  ✓ Binary search for findLastNonOverlapping"
    echo "  ✓ O(n²) → O(n log n) complexity reduction"
    echo "  ✓ Parallel DP for large problems"
    echo "  ✓ Performance comparison with original DP"
    echo ""
    exit 0
else
    echo "❌ TEST FAILED"
    exit 1
fi
