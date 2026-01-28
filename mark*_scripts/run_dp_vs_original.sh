#!/bin/bash
#SBATCH --job-name=dp_vs_orig
#SBATCH --output=/home/users/lz280/dp_vs_original_%j.out
#SBATCH --error=/home/users/lz280/dp_vs_original_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "================================================================================"
echo "DP-MARKStar vs Original MARKStar Performance Comparison"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""
echo "This comprehensive test compares:"
echo "  - Original MARKStar with Greedy correction selection"
echo "  - DP-MARKStar with Phase 1 optimizations (Binary Search + Parallel DP)"
echo ""
echo "Test scales: 5, 7, and 9 flexible residues"
echo ""
echo "Expected runtime: 3-5 hours"
echo "================================================================================"
echo ""

# Use isolated Gradle home
export GRADLE_USER_HOME=/tmp/gradle_comparison_$$
mkdir -p $GRADLE_USER_HOME

# Run comprehensive comparison test
echo "Running comprehensive comparison test..."
echo ""

./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestDPvsOriginal.testComprehensiveComparison" \
    --rerun-tasks 2>&1

echo ""
echo "================================================================================"
echo "Test completed at: $(date)"
echo "================================================================================"

# Cleanup
rm -rf $GRADLE_USER_HOME
