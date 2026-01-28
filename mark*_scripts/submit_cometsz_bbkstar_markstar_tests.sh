#!/bin/bash -l

#SBATCH --job-name=cometsz-bbkstar-markstar
#SBATCH --output=cometsz_bbkstar_markstar_%j.out
#SBATCH --error=cometsz_bbkstar_markstar_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

export GRADLE_USER_HOME=/tmp/$SLURM_JOB_ID

# Print job information
echo "=========================================="
echo "SLURM Job Configuration"
echo "=========================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Job name:      $SLURM_JOB_NAME"
echo "Node:          $SLURM_NODELIST"
echo "CPUs:          $SLURM_CPUS_PER_TASK"
echo "Memory:        $SLURM_MEM_PER_NODE MB ($(($SLURM_MEM_PER_NODE / 1024)) GB)"
echo "Start time:    $(date)"
echo "Working dir:   $(pwd)"
echo "=========================================="
echo ""

# Load Java module (adjust version as needed)
module load Java/17

# Verify Java version
echo "Java version:"
java -version 2>&1
echo ""

# Set JAVA_HOME if needed
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
echo "JAVA_HOME: $JAVA_HOME"
echo ""

# Set JVM memory options to use available memory
# Use 80% of allocated memory for JVM heap
JVM_HEAP_GB=$((SLURM_MEM_PER_NODE * 80 / 100 / 1024))
export GRADLE_OPTS="-Xms${JVM_HEAP_GB}g -Xmx${JVM_HEAP_GB}g -XX:+UseG1GC -XX:MaxGCPauseMillis=200"
echo "JVM Options: $GRADLE_OPTS"
echo ""

# Navigate to project directory
cd /home/users/lz280/IdeaProjects/OSPREY3

# Start resource monitoring in background
(
  echo "=========================================="
  echo "Resource Usage Monitor"
  echo "Test suite: TestCometsZWithBBKStarAndMARKStar"
  echo "=========================================="
  while true; do
    echo "[$(date +%H:%M:%S)] Memory: $(free -h | grep Mem | awk '{print $3 "/" $2}') | CPU Load: $(uptime | awk -F'load average:' '{print $2}')"
    sleep 300  # Report every 5 minutes
  done
) &
MONITOR_PID=$!

# Run CometsZ + BBK* + MARK* performance tests
echo ""
echo "=========================================="
echo "Running CometsZ + BBK* + MARK* Tests"
echo "=========================================="
echo "Test combinations:"
echo "  1. CometsZ + GradientDescent (baseline)"
echo "  2. CometsZ + MARK*"
echo "  3. CometsZ + BBK* + GradientDescent (future)"
echo "  4. CometsZ + BBK* + MARK* (future)"
echo ""
echo "Expected duration: Several hours"
echo "=========================================="
echo ""

# Run tests with info logging
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsZWithBBKStarAndMARKStar.*" --info

# Capture exit code
TEST_EXIT_CODE=$?

# Stop resource monitor
kill $MONITOR_PID 2>/dev/null

# Print completion summary
echo ""
echo "=========================================="
echo "Test Run Summary"
echo "=========================================="
echo "Exit code:     $TEST_EXIT_CODE"
echo "End time:      $(date)"
echo "Node:          $SLURM_NODELIST"
echo ""

# Calculate runtime
if [ -n "$SLURM_JOB_ID" ]; then
  echo "Job statistics will be available via: seff $SLURM_JOB_ID"
fi

# Show final memory usage
echo ""
echo "Final memory usage:"
free -h | grep -E "Mem|Swap"

# Check if tests passed
if [ $TEST_EXIT_CODE -eq 0 ]; then
  echo ""
  echo "✓ All tests PASSED"
  echo ""
  echo "Performance results are in this file."
  echo "Key metrics to extract:"
  echo "  grep 'MARK\* Speedup:' cometsz_bbkstar_markstar_${SLURM_JOB_ID}.out"
  echo "  grep 'Total Time:' cometsz_bbkstar_markstar_${SLURM_JOB_ID}.out"
  echo "  grep 'Execution Time:' cometsz_bbkstar_markstar_${SLURM_JOB_ID}.out"
else
  echo ""
  echo "✗ Tests FAILED with exit code $TEST_EXIT_CODE"
  echo "Check error file: cometsz_bbkstar_markstar_${SLURM_JOB_ID}.err"
fi

echo "=========================================="

# Exit with the test exit code
exit $TEST_EXIT_CODE
