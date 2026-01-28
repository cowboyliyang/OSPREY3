#!/bin/bash
#SBATCH --job-name=osprey_bugfix_test
#SBATCH --output=test_bugfix_%j.out
#SBATCH --error=test_bugfix_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --partition=compsci

# Load Java module
module load Java/17.0.4

# Set working directory
cd /home/users/lz280/IdeaProjects/OSPREY3

echo "========================================"
echo "OSPREY3 SubtreeDOFCache Bug Fix Test"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "========================================"
echo ""

# Clean build
echo "Cleaning previous build..."
./gradlew clean

# Compile the project
echo ""
echo "Compiling project..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo ""
echo "========================================"
echo "Running Bug Fix Test: testComprehensiveComparison"
echo "This test verifies the SubtreeDOFCache fix"
echo "========================================"
echo ""

# Run the specific test that exercises the bug fix
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPvsOriginal.testComprehensiveComparison"

TEST_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Test completed with exit code: $TEST_EXIT_CODE"
echo "End time: $(date)"
echo "========================================"

# Check test results
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✓ TEST PASSED - Bug fix verified!"
    echo ""
    echo "The fix successfully ensures that:"
    echo "  1. Protein ConfSpace uses its own cache (5 positions)"
    echo "  2. Ligand ConfSpace uses its own cache (3 positions)"
    echo "  3. Complex ConfSpace uses its own cache (8 positions)"
    echo "  4. No position index incompatibility errors occur"
else
    echo ""
    echo "✗ TEST FAILED - Check error logs"
    echo ""
fi

# Copy test results to analysis directory
echo ""
echo "Copying test output to analysis directory..."
cp build/reports/tests/test/classes/edu.duke.cs.osprey.markstar.TestDPvsOriginal.html \
   /home/users/lz280/"260114 triple correction analysis"/test_results_bugfix_${SLURM_JOB_ID}.html 2>/dev/null || echo "HTML report not found"

# Look for test XML output
if [ -f "build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml" ]; then
    cp build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml \
       /home/users/lz280/"260114 triple correction analysis"/test_results_bugfix_${SLURM_JOB_ID}.xml
    echo "Test XML results copied"
fi

echo ""
echo "========================================"
echo "Job complete!"
echo "========================================"

exit $TEST_EXIT_CODE
