#!/bin/bash
#SBATCH --job-name=comets-z-full
#SBATCH --output=comets_full_%j.out
#SBATCH --error=comets_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""

# Load Java
module load Java/17

# Set JAVA_HOME
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))

# Navigate to project
cd /home/users/lz280/IdeaProjects/OSPREY3

# Run full 2RL0 test (8 positions, 25 sequences - this will take a while!)
echo "Running full 2RL0 COMETS-Z performance test (this may take hours)..."
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsVsCometsZPerformance.compare2RL0Full"

echo "End time: $(date)"
