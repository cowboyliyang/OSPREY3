#!/bin/bash
#SBATCH --job-name=mskstar-tests
#SBATCH --output=mskstar_%j.out
#SBATCH --error=mskstar_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
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

# Run all original MSKStar tests
echo "Running all MSKStar tests..."
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMSKStar.*"

echo "End time: $(date)"
