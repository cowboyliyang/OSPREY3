#!/bin/bash
#SBATCH --job-name=triple_cache_test
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/test_output_%j.log
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/test_error_%j.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=compsci

# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Starting time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Change to project directory
cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean all locks first
echo "=== Cleaning Gradle locks ==="
find ~/.gradle -name "*.lock" -delete 2>/dev/null
find .gradle -name "*.lock" -delete 2>/dev/null
find buildSrc/.gradle -name "*.lock" -delete 2>/dev/null

# Force clean build to ensure new code is compiled
echo ""
echo "=== Force cleaning build cache ==="
rm -rf build .gradle/*/classes buildSrc/.gradle buildSrc/build
echo "Build cache cleared"

# Stop any existing gradle daemons
echo ""
echo "=== Stopping existing gradle daemons ==="
./gradlew --stop 2>/dev/null || true
sleep 2
echo ""

# Run the test
echo "Running triple DOF cache comparison test..."
echo "This will run 3 versions: Original, Subtree Cache, Subtree+Triple Cache"
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestTripleDOFCache.testMinimizationTiming6Mutable3Flexible" --no-daemon --info

echo ""
echo "Test completed at: $(date)"
echo ""
echo "Results saved to:"
echo "  Test XML: build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestTripleDOFCache.xml"
echo "  Test HTML: build/reports/tests/test/classes/edu.duke.cs.osprey.markstar.TestTripleDOFCache.html"
