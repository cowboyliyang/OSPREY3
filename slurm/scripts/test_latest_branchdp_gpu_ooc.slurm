#!/bin/bash
# Two-GPU regression for the current branch-DP OOC planner and lambda-box
# dispatch.  All compile/test execution stays inside this Slurm job.
#SBATCH --job-name=packstar_gpu_ooc
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/gpu_ooc_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/gpu_ooc_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:2

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/gpu-ooc-build-${SLURM_JOB_ID}
GRADLE_ROOT=/usr/xtmp/lz280/packstar_migration/gradle-home

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT"
export GRADLE_USER_HOME="$GRADLE_ROOT"

if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load cuda/12.8 >/dev/null 2>&1 \
        || module load cuda/12.4 >/dev/null 2>&1 \
        || module load cuda >/dev/null 2>&1 \
        || true
fi
# A conda CUDA environment may export an older nvcc through NVCC/CUDA_HOME;
# the Gradle task has fixed system-toolkit fallbacks that support the SM90
# fatbin target, so do not let that environment override them.
unset NVCC CUDA_HOME

cd "$REPO"
mkdir -p /usr/xtmp/lz280/packstar_migration
nvidia-smi --query-gpu=index,name,memory.total,memory.free \
    --format=csv,noheader

./gradlew --no-daemon --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    compileCuda_dp compileTestJava

./gradlew --no-daemon --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx4g -DtestMaxHeap=8g \
    test \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.outOfCorePlannerAndPackingBoundBothRowAndLambdaAxes_withoutGpu \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.outOfCoreStablePartialMergeHandlesExtremeNanAndNegativeInfinity_withoutGpu \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCoreExtremeNanNegativeInfinityAndAllocationAudit \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_bilateralRowAndLambdaTilingMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_fileBackedChildrenAndOutputMatchJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuOutOfCore_multiGpuMatchesJava

echo "branch-DP GPU/OOC regression passed job=$SLURM_JOB_ID"
