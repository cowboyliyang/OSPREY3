#!/bin/bash
#SBATCH --job-name=osprey_perf_test
#SBATCH --output=test_perf_%j.out
#SBATCH --error=test_perf_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --partition=compsci

# Load Java module
module load Java/17.0.4

# Set working directory
cd /home/users/lz280/IdeaProjects/OSPREY3

echo "========================================"
echo "OSPREY3 Performance Comparison Test"
echo "Testing: Original vs Phase 1+2 (with bugfix)"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "========================================"
echo ""

echo "Running performance comparison test..."
echo "This will compare Original vs Phase 1+2 with the bugfix applied"
echo ""

# Run the comprehensive performance test
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPvsOriginal.testAllPhasesIntegrated"

TEST_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Test completed with exit code: $TEST_EXIT_CODE"
echo "End time: $(date)"
echo "========================================"

# Copy test results to analysis directory
echo ""
echo "Copying test output to analysis directory..."
if [ -f "build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml" ]; then
    cp build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml \
       /home/users/lz280/"260114 triple correction analysis"/test_performance_${SLURM_JOB_ID}.xml
    echo "Test XML results copied"
fi

# Extract performance summary from output
echo ""
echo "========================================"
echo "Extracting Performance Summary..."
echo "========================================"
grep -A 20 "FINAL PERFORMANCE SUMMARY" build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml 2>/dev/null || echo "Summary not found in XML"

echo ""
echo "========================================"
echo "Job complete!"
echo "========================================"

exit $TEST_EXIT_CODE
