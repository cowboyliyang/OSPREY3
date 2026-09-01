#!/bin/bash
# Targeted validation for the latest migrated PackStar and exact higher-order
# DP/GPU paths.  Run through Slurm so all computation stays off the login node.
#SBATCH --job-name=packstar_migration_test
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/test_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/test_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/gradle-build
GRADLE_ROOT=/usr/xtmp/lz280/packstar_migration/gradle-home

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT"
export GRADLE_USER_HOME="$GRADLE_ROOT"

command -v module >/dev/null 2>&1 && {
    module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
}
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

cd "$REPO"
./gradlew test \
    --no-daemon \
    --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/packstar_test_filter.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g \
    -DtestMaxHeap=16g \
    --tests edu.duke.cs.osprey.packstar.TestPackStarConfig \
    --tests edu.duke.cs.osprey.packstar.TestPackStarLogSpaceBounds \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleEtaCorrections \
    --tests edu.duke.cs.osprey.branchdp.TestHigherOrderFullDP \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_nonLeafRegularMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_childSlicedMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_hybridMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_outOfCoreMatchesJava
