#!/bin/bash
#SBATCH --job-name=dp_complete
#SBATCH --output=/home/users/lz280/dp_complete_%j.out
#SBATCH --error=/home/users/lz280/dp_complete_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "================================================================================"
echo "COMPLETE DP-MARKSTAR EVALUATION: All Phases Integrated"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""
echo "This comprehensive test compares:"
echo "  1. Original MARKStar (Greedy correction selection)"
echo "  2. Phase 1: DP with Binary Search + Parallel DP"
echo "  3. Phase 1+2: DP + Subtree DOF Cache with Branch Decomposition"
echo ""
echo "Test scales: 7 and 9 flexible residues"
echo ""
echo "Expected runtime: 5-7 hours"
echo "================================================================================"
echo ""

# Use isolated Gradle home
export GRADLE_USER_HOME=/tmp/gradle_complete_$$
mkdir -p $GRADLE_USER_HOME

# Run complete integrated test
echo "Running Phase 1+2 integrated test..."
echo ""

./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestDPvsOriginal.testAllPhasesIntegrated" \
    --rerun-tasks 2>&1

echo ""
echo "================================================================================"
echo "Test completed at: $(date)"
echo "================================================================================"

# Cleanup
rm -rf $GRADLE_USER_HOME
