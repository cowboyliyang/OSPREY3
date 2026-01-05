#!/bin/bash
#SBATCH --job-name=comets-z-small
#SBATCH --output=comets_small_%j.out
#SBATCH --error=comets_small_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
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

# Run small tests only
echo "Running small COMETS-Z performance tests..."
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsVsCometsZPerformance.compareSmall2RL0" \
               --tests "edu.duke.cs.osprey.kstar.TestCometsVsCometsZPerformance.compareSmall2RL0BoundedMemory"

echo "End time: $(date)"
