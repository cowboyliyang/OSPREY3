#!/bin/bash
#SBATCH --job-name=dp_phase1_opt
#SBATCH --output=phase1_opt_%j.out
#SBATCH --error=phase1_opt_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compsci

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=== DP-MARKStar Phase 1 Optimizations Test ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Starting at: $(date)"
echo ""

# Run the test directly
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestDPOptimizations.testPhase1Optimizations" 2>&1

echo ""
echo "Finished at: $(date)"
