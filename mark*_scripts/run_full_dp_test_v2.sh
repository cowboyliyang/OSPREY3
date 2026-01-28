#!/bin/bash
#SBATCH --job-name=dp_full_v2
#SBATCH --output=dp_full_v2_%j.log
#SBATCH --error=dp_full_v2_%j.err
#SBATCH --time=30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP-MARKStar Complete Feature Test V2"
echo "All enhancements + Multiple scales"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 16GB, CPUs: 4, Time limit: 30min"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Stop any existing gradle daemons
./gradlew --stop || true
sleep 3

echo "Compiling with all DP enhancements..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful"
echo ""

echo "=========================================="
echo "Test Suite:"
echo "  1. Small scale (2 residues) - verify features"
echo "  2. Medium scale (5 residues) - verify performance"
echo "=========================================="
echo ""

# Test 1: Small scale feature verification
echo "==================== TEST 1: Small Scale ===================="
echo "Verifying all DP enhancements (2 flexible residues)..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"

TEST1_RESULT=$?

if [ $TEST1_RESULT -eq 0 ]; then
    echo ""
    echo "✓ TEST 1 PASSED: Feature verification successful"
else
    echo ""
    echo "❌ TEST 1 FAILED: Feature verification failed"
fi

echo ""
echo "==================== TEST 2: Medium Scale ===================="
echo "Testing DP performance (5 flexible residues)..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPComplete.testMediumScale"

TEST2_RESULT=$?

if [ $TEST2_RESULT -eq 0 ]; then
    echo ""
    echo "✓ TEST 2 PASSED: Medium scale test successful"
else
    echo ""
    echo "❌ TEST 2 FAILED: Medium scale test failed"
fi

echo ""
echo "==================== TEST 3: Complete Suite ===================="
echo "Running complete feature test..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPComplete.testCompleteFeatures"

TEST3_RESULT=$?

if [ $TEST3_RESULT -eq 0 ]; then
    echo ""
    echo "✓ TEST 3 PASSED: Complete suite successful"
else
    echo ""
    echo "❌ TEST 3 FAILED"
fi

echo ""
echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo "End time: $(date)"
echo ""

TOTAL_PASSED=0
TOTAL_FAILED=0

if [ $TEST1_RESULT -eq 0 ]; then TOTAL_PASSED=$((TOTAL_PASSED+1)); else TOTAL_FAILED=$((TOTAL_FAILED+1)); fi
if [ $TEST2_RESULT -eq 0 ]; then TOTAL_PASSED=$((TOTAL_PASSED+1)); else TOTAL_FAILED=$((TOTAL_FAILED+1)); fi
if [ $TEST3_RESULT -eq 0 ]; then TOTAL_PASSED=$((TOTAL_PASSED+1)); else TOTAL_FAILED=$((TOTAL_FAILED+1)); fi

echo "Tests Passed: $TOTAL_PASSED / 3"
echo "Tests Failed: $TOTAL_FAILED / 3"
echo ""

if [ $TOTAL_FAILED -eq 0 ]; then
    echo "✓✓✓ ALL TESTS PASSED ✓✓✓"
    echo ""
    echo "Verified Features:"
    echo "  ✓ Backtracking (getSelectedCorrectionsFromDP)"
    echo "  ✓ Detailed statistics tracking"
    echo "  ✓ Enhanced printDPStats() output"
    echo "  ✓ DP vs Greedy comparison"
    echo "  ✓ Multi-scale performance"
    echo ""
    exit 0
else
    echo "❌ $TOTAL_FAILED TEST(S) FAILED"
    exit 1
fi
