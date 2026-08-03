#!/bin/bash
#SBATCH --job-name=branchdp_gpu_tests
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_tests_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/branchdp_gpu_tests_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

echo "=============================================="
echo "Branch-DP targeted GPU tests"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "Start: $(date)"
echo "Java: $(java -version 2>&1 | head -n 1)"
echo "=============================================="

./gradlew test \
    --no-daemon \
    --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g \
    -DtestMaxHeap=16g \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP \
    --tests edu.duke.cs.osprey.branchdp.TestComputeFullDPSharded \
    --tests edu.duke.cs.osprey.branchdp.TestDPTable

echo ""
echo "End: $(date)"
