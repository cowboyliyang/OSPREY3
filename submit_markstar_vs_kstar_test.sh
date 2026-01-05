#!/bin/bash
#SBATCH --job-name=markstar_vs_kstar_test
#SBATCH --output=markstar_vs_kstar_%j.out
#SBATCH --error=markstar_vs_kstar_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compsci

# Load Java module if needed
# module load java/11

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=========================================="
echo "Testing MARK* vs K* Partition Function"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run the test
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction.testSmallSystem3Flex" --info

echo ""
echo "=========================================="
echo "End time: $(date)"
echo "=========================================="

# Print memory usage
echo ""
echo "==========================================
Test Run Summary
==========================================
Exit code:     $?
End time:      $(date)
Node:          $(hostname)

Job statistics will be available via: seff $SLURM_JOB_ID

Final memory usage:"
free -h
echo "=========================================="
