#!/bin/bash
#SBATCH --job-name=dp_mark_large
#SBATCH --output=dp_comparison_%j.log
#SBATCH --error=dp_comparison_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

echo "=========================================="
echo "DP-MARKStar vs Greedy Comparison on 2RL0"
echo "LARGE SCALE TEST (9 flexible residues)"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Memory: 32GB, CPUs: 4, Time limit: 6h"
echo ""

cd /home/users/lz280/IdeaProjects/OSPREY3

# Stop any existing gradle daemons to avoid lock conflicts
./gradlew --stop || true
sleep 5

# Run the comparison test
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDP_MARKStar.test2RL0_GreedyVsDP"

echo ""
echo "End time: $(date)"
echo "Test completed!"
