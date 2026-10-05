#!/bin/bash
#SBATCH --job-name=ps-performance-check
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --constraint=a5000
#SBATCH --gres=gpu:a5000:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_performance_20261003/check_${SLURM_JOB_ID}
# Source/fixture snapshot ~130 MiB / several thousand files; build and reports <2 GiB expected.
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
export PACKSTAR_CLASSPATH_FILE="$RUN_ROOT/test_classpath.txt"
cd "$REPO"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib examples buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
git diff --binary HEAD > "$RUN_ROOT/source.patch"
git ls-files --others --exclude-standard > "$RUN_ROOT/untracked_sources.txt"
cd "$RUN_ROOT/source"
./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    -I slurm/scripts/gradle_xtmp.init.gradle -I slurm/scripts/verify_ccd.init.gradle \
    -I slurm/h200/classpath.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" -DtestMaxHeap=16g \
    --tests 'edu.duke.cs.osprey.packstar.TestPackStar*' \
    --tests edu.duke.cs.osprey.tools.TestByteBoundedCache \
    --tests edu.duke.cs.osprey.parallelism.TestTaskCapacity \
    --tests edu.duke.cs.osprey.branchdp.TestChildFold \
    --tests edu.duke.cs.osprey.branchdp.TestComputeFullDPSharded \
    --tests edu.duke.cs.osprey.branchdp.TestHigherOrderFullDP \
    --tests edu.duke.cs.osprey.branchdp.TestRootedTreeMemoryRelease \
    --tests edu.duke.cs.osprey.branchdp.TestGpuPhase1Sampling \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.samplingPreparationIsBoundedAndScopedToAnExplicitProposalEpoch \
    writeExperimentClasspath
echo "PERFORMANCE_VALIDATION=$RUN_ROOT"
