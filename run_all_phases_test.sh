#!/bin/bash
#SBATCH --job-name=dp_all_phases
#SBATCH --output=/home/users/lz280/all_phases_%j.out
#SBATCH --error=/home/users/lz280/all_phases_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=== DP-MARKStar All Phases Test ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Starting at: $(date)"
echo ""
echo "This test will run:"
echo "  - testPhase1Optimizations: Phase 1 binary search + parallel DP"
echo "  - testPhase2BranchDecomposition: Phase 2 branch decomposition"
echo "  - testPhase3DPTrie: Phase 3 DP-Trie integration"
echo "  - testAllPhasesCombined: Comprehensive test"
echo ""

# Use isolated Gradle home
export GRADLE_USER_HOME=/tmp/gradle_all_phases_$$
mkdir -p $GRADLE_USER_HOME

# Run all tests
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPOptimizations" --rerun-tasks 2>&1

# Cleanup
rm -rf $GRADLE_USER_HOME

echo ""
echo "Finished at: $(date)"
