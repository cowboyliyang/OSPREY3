#!/bin/bash
#SBATCH --job-name=dp_verify_v2
#SBATCH --output=dp_verify_v2_%j.log
#SBATCH --error=dp_verify_v2_%j.err
#SBATCH --time=30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP Enhancements Verification V2"
echo "Using isolated Gradle user home"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Use isolated Gradle home to avoid conflicts
export GRADLE_USER_HOME=/tmp/gradle_verify_$$
mkdir -p $GRADLE_USER_HOME

echo "Using isolated Gradle home: $GRADLE_USER_HOME"
echo ""

echo "Compiling with DP enhancements..."
./gradlew compileJava compileTestJava

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    rm -rf $GRADLE_USER_HOME
    exit 1
fi

echo "✓ Compilation successful"
echo ""

# Run small-scale verification test
echo "Running DP enhancements verification test..."
echo ""

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPEnhancements.testDPEnhancementsSmall"

TEST_RESULT=$?

# Cleanup
rm -rf $GRADLE_USER_HOME

echo ""
echo "End time: $(date)"

if [ $TEST_RESULT -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓✓✓ VERIFICATION SUCCESSFUL ✓✓✓"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "❌ VERIFICATION FAILED"
    echo "=========================================="
    exit 1
fi
