#!/bin/bash
# Build the current dirty working tree once for the rescue array.  All Gradle
# outputs and caches are redirected to xtmp; array workers only read this tree.
#SBATCH --job-name=ps_rescue_build
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=${BUILD_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/current_build}
GRADLE_CACHE=${GRADLE_CACHE:-/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache}
INIT="$REPO/slurm/scripts/gradle_xtmp.init.gradle"
CP_INIT="$REPO/slurm/scripts/print_test_classpath.init.gradle"
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}

mkdir -p "$BUILD_ROOT" "$GRADLE_CACHE"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT/gradle_build"
export GRADLE_USER_HOME="$GRADLE_CACHE"

cd "$REPO"
echo "build_job=$SLURM_JOB_ID start=$(date --iso-8601=seconds)"
echo "repo=$REPO"
echo "build_root=$BUILD_ROOT"
git status --short
git rev-parse HEAD

./gradlew -I "$INIT" testClasses --no-daemon --max-workers="$SLURM_CPUS_PER_TASK"
./gradlew -q -I "$INIT" -I "$CP_INIT" \
    printTestClasspath --no-daemon --max-workers=1 \
    > "$BUILD_ROOT/test_classpath.txt"

sha256sum \
    src/main/java/edu/duke/cs/osprey/packstar/PackStarEstimator.java \
    src/main/java/edu/duke/cs/osprey/packstar/PackStarTripleEtaCorrections.java \
    src/main/java/edu/duke/cs/osprey/packstar/PackStarProposalLearning.java \
    src/test/java/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.java \
    "$BUILD_ROOT/test_classpath.txt" \
    > "$BUILD_ROOT/source_and_classpath.sha256"
{
    printf 'buildJob\t%s\n' "$SLURM_JOB_ID"
    printf 'node\t%s\n' "$SLURMD_NODENAME"
    printf 'gitHead\t%s\n' "$(git rev-parse HEAD)"
    printf 'buildRoot\t%s\n' "$BUILD_ROOT"
    printf 'startedAt\t%s\n' "$(date --iso-8601=seconds)"
} > "$BUILD_ROOT/build_manifest.tsv"

echo "build_job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
