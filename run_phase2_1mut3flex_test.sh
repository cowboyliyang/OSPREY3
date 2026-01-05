#!/bin/bash
#SBATCH -p compsci
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 20
#SBATCH -t 30:00
#SBATCH --mail-type=all
#SBATCH --mail-user=lz280@duke.edu

echo "========================================="
echo "Phase 2 Test: 1 Mutable + 3 Flexible"
echo "Total: 4 positions (G648 mutable + 3 wild-type)"
echo "Original vs Phase 1 vs Phase 1+2"
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "========================================="

module load java/17

cd /home/users/lz280/IdeaProjects/OSPREY3

echo ""
echo "Running comprehensive comparison on 3 flexible residues..."
echo "Configuration: buildConfSpace(3) → 1 mutable (G648) + 3 flexible"
echo "This will test:"
echo "  1. Original Greedy (baseline)"
echo "  2. Phase 1+2 (DP + TRUE Subtree Caching)"
echo "  (Phase 1 skipped)"
echo ""

./gradlew test --tests TestDPvsOriginal.testComprehensiveComparison 2>&1

echo ""
echo "========================================="
echo "Test completed with exit code: $?"
echo "Date: $(date)"
echo "========================================="
