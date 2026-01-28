#!/bin/bash
#SBATCH --job-name=dp_phase1_test
#SBATCH --output=/home/users/lz280/phase1_test_%j.out
#SBATCH --error=/home/users/lz280/phase1_test_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=== DP-MARKStar Phase 1 Optimizations Test ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Starting at: $(date)"
echo ""

# Use isolated Gradle home
export GRADLE_USER_HOME=/tmp/gradle_test_$$
mkdir -p $GRADLE_USER_HOME

# Clean and run tests
./gradlew clean
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPOptimizations.testPhase1Optimizations" --rerun-tasks 2>&1

# Cleanup
rm -rf $GRADLE_USER_HOME

echo ""
echo "Finished at: $(date)"
