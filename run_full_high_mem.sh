#!/bin/bash
#SBATCH --job-name=dp_highmem
#SBATCH --output=dp_highmem_%j.log
#SBATCH --error=dp_highmem_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compsci

echo "=========================================="
echo "DP-MARKStar High Memory Complete Test"
echo "Maximum memory for optimal DP performance"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 64GB, CPUs: 8, Time limit: 1h"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Stop any existing gradle daemons
./gradlew --stop || true
sleep 3

echo "Compiling with all DP enhancements..."
echo "(Using high memory for DP caching and optimization)"
echo ""

./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful"
echo ""

echo "=========================================="
echo "High-Memory Test Suite"
echo "  Memory strategy: Space-for-Time tradeoff"
echo "  DP benefits from more memory for:"
echo "    - Larger DP arrays"
echo "    - More backtracking history"
echo "    - Better cache utilization"
echo "=========================================="
echo ""

# Test 1: Small scale
echo "==================== TEST 1: Small (2 res) ===================="
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"
TEST1=$?
echo "TEST 1 result: $TEST1"
echo ""

# Test 2: Medium scale
echo "==================== TEST 2: Medium (5 res) ===================="
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPComplete.testMediumScale"
TEST2=$?
echo "TEST 2 result: $TEST2"
echo ""

# Test 3: Complete features
echo "==================== TEST 3: Complete Features ===================="
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPComplete.testCompleteFeatures"
TEST3=$?
echo "TEST 3 result: $TEST3"
echo ""

# Summary
echo ""
echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo "End time: $(date)"
echo ""

PASSED=0
FAILED=0

if [ $TEST1 -eq 0 ]; then
    echo "✓ TEST 1 PASSED (Small scale)"
    PASSED=$((PASSED+1))
else
    echo "❌ TEST 1 FAILED"
    FAILED=$((FAILED+1))
fi

if [ $TEST2 -eq 0 ]; then
    echo "✓ TEST 2 PASSED (Medium scale)"
    PASSED=$((PASSED+1))
else
    echo "❌ TEST 2 FAILED"
    FAILED=$((FAILED+1))
fi

if [ $TEST3 -eq 0 ]; then
    echo "✓ TEST 3 PASSED (Complete features)"
    PASSED=$((PASSED+1))
else
    echo "❌ TEST 3 FAILED"
    FAILED=$((FAILED+1))
fi

echo ""
echo "Results: $PASSED passed, $FAILED failed out of 3 tests"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓✓✓ ALL TESTS PASSED WITH HIGH MEMORY ✓✓✓"
    echo ""
    echo "Verified DP Enhancements:"
    echo "  ✓ Backtracking (getSelectedCorrectionsFromDP)"
    echo "  ✓ Detailed statistics (selection rates)"
    echo "  ✓ Enhanced printDPStats() output"
    echo "  ✓ Multi-scale testing (2, 5 residues)"
    echo ""
    echo "Memory Configuration:"
    echo "  Allocated: 64GB"
    echo "  CPUs: 8"
    echo "  DP space-for-time optimization: ENABLED"
    echo ""
    exit 0
else
    echo "❌ SOME TESTS FAILED"
    echo "Check logs for details"
    exit 1
fi
