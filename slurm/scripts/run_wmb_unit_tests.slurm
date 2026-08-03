#!/bin/bash
#SBATCH --job-name=wmb_unit
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/wmb_unit_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/wmb_unit_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

echo "=============================================="
echo "WMB unit tests"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "Java: $(java -version 2>&1 | head -n 1)"
echo "=============================================="

./gradlew compileJava compileTestJava test --tests edu.duke.cs.osprey.wmb.TestWeightedMiniBucket --no-daemon

echo ""
echo "End: $(date)"
