#!/bin/bash
#SBATCH --job-name=branchdp_gpu_hybrid
#SBATCH --output=/usr/xtmp/lz280/bench_comparison/logs/branchdp_gpu_hybrid_%j.out
#SBATCH --error=/usr/xtmp/lz280/bench_comparison/logs/branchdp_gpu_hybrid_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:2

set -euo pipefail

cd /usr/xtmp/lz280/osprey3_test_validate

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
module load cuda/cuda-12.1

echo "node=$SLURM_NODELIST gpu=${CUDA_VISIBLE_DEVICES:-unset} start=$(date)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader

./gradlew test --offline --no-daemon --max-workers=1 \
    -Dorg.gradle.jvmargs=-Xmx2g -DtestMaxHeap=4g \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuHybridResidentStreamed_matchesJava_twoChildren \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuHybridResidentStreamed_multiGpuMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.hybridPlannerGatherAndMEnumeration_withoutGpu \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.hybridPlannerChooses3buaResidentStreamedShape_withoutAllocatingTables \
    --tests edu.duke.cs.osprey.branchdp.TestDPTable

echo "end=$(date)"
