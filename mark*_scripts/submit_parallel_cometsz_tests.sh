#!/bin/bash -l

#SBATCH --job-name=parallel-cometsz
#SBATCH --output=parallel_cometsz_%j.out
#SBATCH --error=parallel_cometsz_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

export GRADLE_USER_HOME=/tmp/$SLURM_JOB_ID

# Print job information
echo "=========================================="
echo "Parallel COMETS-Z Tests Configuration"
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

# Load Java module
module load Java/17

# Verify Java version
echo "Java version:"
java -version 2>&1
echo ""

# Set JAVA_HOME
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
echo "JAVA_HOME: $JAVA_HOME"
echo ""

# Set JVM memory options
# Allocate 40GB per test (for 2 parallel tests = 80GB total)
export GRADLE_OPTS="-Xmx40g -Xms40g -XX:+UseG1GC -XX:MaxGCPauseMillis=200"
echo "JVM Options: $GRADLE_OPTS"
echo ""

# Navigate to project directory
cd /home/users/lz280/IdeaProjects/OSPREY3

# Start resource monitoring in background
(
  echo "=========================================="
  echo "Resource Usage Monitor"
  echo "=========================================="
  while true; do
    timestamp=$(date +%H:%M:%S)
    mem_usage=$(free -h | grep Mem | awk '{print $3 "/" $2}')
    cpu_load=$(uptime | awk -F'load average:' '{print $2}')
    echo "[$timestamp] Memory: $mem_usage | CPU Load: $cpu_load"
    sleep 300  # Report every 5 minutes
  done
) &
MONITOR_PID=$!

# Build project first (avoid parallel build conflicts)
echo "=========================================="
echo "Building project (one-time build)..."
echo "=========================================="
./gradlew build -x test --info
BUILD_EXIT_CODE=$?

if [ $BUILD_EXIT_CODE -ne 0 ]; then
  echo ""
  echo "✗ Build FAILED with exit code $BUILD_EXIT_CODE"
  kill $MONITOR_PID 2>/dev/null
  exit $BUILD_EXIT_CODE
fi

echo ""
echo "✓ Build completed successfully"
echo ""

# Run tests in parallel
echo "=========================================="
echo "Running tests in parallel..."
echo "=========================================="
echo "Tests to run:"
echo "  1. TestCometsZWithMARKStarPerformance"
echo "  2. TestCometsZWithBBKStarAndMARKStar"
echo "=========================================="
echo ""

# Create output directory for individual test logs
mkdir -p test_logs_${SLURM_JOB_ID}

# Run tests in background
echo "[$(date +%H:%M:%S)] Starting TestCometsZWithMARKStarPerformance..."
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsZWithMARKStarPerformance.*" --info \
  > test_logs_${SLURM_JOB_ID}/cometsz_markstar.log 2>&1 &
PID1=$!

echo "[$(date +%H:%M:%S)] Starting TestCometsZWithBBKStarAndMARKStar..."
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsZWithBBKStarAndMARKStar.*" --info \
  > test_logs_${SLURM_JOB_ID}/cometsz_bbkstar_markstar.log 2>&1 &
PID2=$!

echo ""
echo "All tests started. PIDs:"
echo "  TestCometsZWithMARKStarPerformance: $PID1"
echo "  TestCometsZWithBBKStarAndMARKStar:  $PID2"
echo ""
echo "Waiting for all tests to complete..."
echo "This may take several hours..."
echo ""

# Wait for all tests and track completion
wait $PID1
EXIT1=$?
echo "[$(date +%H:%M:%S)] TestCometsZWithMARKStarPerformance completed with exit code: $EXIT1"

wait $PID2
EXIT2=$?
echo "[$(date +%H:%M:%S)] TestCometsZWithBBKStarAndMARKStar completed with exit code: $EXIT2"

# Stop resource monitor
kill $MONITOR_PID 2>/dev/null

# Print test results summary
echo ""
echo "=========================================="
echo "Test Results Summary"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Results:"
echo "  1. TestCometsZWithMARKStarPerformance: $([ $EXIT1 -eq 0 ] && echo '✓ PASSED' || echo '✗ FAILED (exit code '$EXIT1')')"
echo "  2. TestCometsZWithBBKStarAndMARKStar:  $([ $EXIT2 -eq 0 ] && echo '✓ PASSED' || echo '✗ FAILED (exit code '$EXIT2')')"
echo ""

# Show log file locations
echo "Individual test logs saved in:"
echo "  test_logs_${SLURM_JOB_ID}/cometsz_markstar.log"
echo "  test_logs_${SLURM_JOB_ID}/cometsz_bbkstar_markstar.log"
echo ""

# Calculate runtime
if [ -n "$SLURM_JOB_ID" ]; then
  echo "Job statistics will be available via: seff $SLURM_JOB_ID"
fi

# Show final memory usage
echo ""
echo "Final memory usage:"
free -h | grep -E "Mem|Swap"
echo ""

# Overall status
if [ $EXIT1 -eq 0 ] && [ $EXIT2 -eq 0 ]; then
  echo "✓ ALL TESTS PASSED"
  echo "=========================================="
  exit 0
else
  echo "✗ SOME TESTS FAILED"
  echo "Check individual log files for details."
  echo "=========================================="
  exit 1
fi
