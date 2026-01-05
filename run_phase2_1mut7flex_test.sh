#!/bin/bash
#SBATCH -p compsci
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=lz280@duke.edu

=========================================
Phase 2 Test: 1 Mutable + 7 Flexible
Original vs Phase 1 vs Phase 1+2
Date: $(date)
Job ID: $SLURM_JOB_ID
Node: $(hostname)
=========================================

module load java/17

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "Running comprehensive comparison on 7 and 9 residues (1 mutable + flexible)..."
echo "This will test:"
echo "  1. Original Greedy (baseline)"
echo "  2. Phase 1 (DP optimizations)"
echo "  3. Phase 1+2 (DP + TRUE Subtree Caching)"
echo ""

./gradlew test --tests TestDPvsOriginal.testComprehensiveComparison 2>&1

echo ""
echo "=========================================
Test completed with exit code: $?
Date: $(date)
========================================="
