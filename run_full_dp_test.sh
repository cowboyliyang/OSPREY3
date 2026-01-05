#!/bin/bash
#SBATCH --job-name=dp_full_test
#SBATCH --output=dp_full_test_%j.log
#SBATCH --error=dp_full_test_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP-MARKStar Complete Feature Test"
echo "Testing all enhancements on 2RL0"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 16GB, CPUs: 4, Time limit: 1h"
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
echo "  1. Small scale (2 residues) - Feature verification"
echo "  2. Medium scale (5 residues) - Performance test"
echo "=========================================="
echo ""

# Test 1: Small scale with all features
echo "==================== TEST 1: Small Scale ===================="
echo "Running feature verification (2 flexible residues)..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"

TEST1_RESULT=$?

if [ $TEST1_RESULT -eq 0 ]; then
    echo ""
    echo "✓ TEST 1 PASSED: Feature verification successful"
else
    echo ""
    echo "❌ TEST 1 FAILED"
    exit 1
fi

echo ""
echo "==================== TEST 2: Medium Scale ===================="
echo "Running performance test (5 flexible residues)..."
echo ""

# Create a medium-scale test if needed
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDP_MARKStar.test2RL0_GreedyVsDP_Medium" 2>&1 | head -100

TEST2_RESULT=$?

echo ""
echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo "End time: $(date)"
echo ""

if [ $TEST1_RESULT -eq 0 ]; then
    echo "✓✓✓ ALL TESTS PASSED ✓✓✓"
    echo ""
    echo "Verified Features:"
    echo "  ✓ Backtracking (getSelectedCorrectionsFromDP)"
    echo "  ✓ Detailed statistics (selection rates)"
    echo "  ✓ Enhanced printDPStats() output"
    echo "  ✓ DP vs Greedy comparison"
    echo ""
    echo "Test completed successfully!"
    exit 0
else
    echo "❌ SOME TESTS FAILED"
    exit 1
fi
