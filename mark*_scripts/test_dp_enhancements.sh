#!/bin/bash
# Quick local test of DP enhancements on small system
# Run this after SLURM job completes to verify new features

echo "=========================================="
echo "Testing DP Enhancements (Local)"
echo "Small scale: 2 flexible residues"
echo "=========================================="
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Wait for any running gradle processes to finish
echo "Waiting for SLURM job to release Gradle lock..."
while pgrep -f "gradle.*TestDP_MARKStar" > /dev/null; do
    echo "  Gradle still running (SLURM job)... waiting 30s"
    sleep 30
done

echo "✓ Gradle lock released"
echo ""

# Stop any daemons
echo "Stopping Gradle daemons..."
./gradlew --stop || true
sleep 3

# Compile with new changes
echo ""
echo "Compiling with DP enhancements..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful"
echo ""

# Run quick test
echo "Running small-scale test to verify new features..."
echo ""

timeout 180 ./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓ DP Enhancements Test PASSED"
    echo "=========================================="
    echo ""
    echo "Verified:"
    echo "  ✓ Backtracking function (getSelectedCorrectionsFromDP)"
    echo "  ✓ Statistics tracking (corrections considered/selected)"
    echo "  ✓ Enhanced printDPStats() formatting"
    echo ""
else
    echo ""
    echo "=========================================="
    echo "❌ Test failed or timed out"
    echo "=========================================="
    echo "Check output above for errors"
    exit 1
fi
