#!/bin/bash
# Snapshot source before building so caches, generated kernels and all test
# artifacts remain on xtmp, including Gradle's buildSrc outputs.
#SBATCH --job-name=packstar_forward_test
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=grisman
#SBATCH --account=grisman
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_forward_20260913/$SLURM_JOB_ID
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export PACKSTAR_TEST_CLASSPATH_FILE="$RUN_ROOT/test_classpath.txt"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
cd "$REPO"
# These source directories total about 60 MiB before excluding build caches.
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
# Keep the build-logic package explicit.  On the previous snapshot this
# directory was omitted even though the rest of buildSrc was copied, which
# made docs.kt fail with unresolved osprey.build references before tests ran.
mkdir -p "$RUN_ROOT/source/buildSrc/src/main/kotlin/osprey/build"
rsync -a buildSrc/src/main/kotlin/osprey/build/ \
    "$RUN_ROOT/source/buildSrc/src/main/kotlin/osprey/build/"
test -f "$RUN_ROOT/source/buildSrc/src/main/kotlin/osprey/build/common.kt"
cd "$RUN_ROOT/source"
./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    --init-script slurm/scripts/forward_classpath.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    -DtestMaxHeap=8g \
    --tests edu.duke.cs.osprey.packstar.TestPackStarForwardTripleSelection \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleDecompositionCosts \
    --tests edu.duke.cs.osprey.packstar.TestPackStarProposalLearning \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleEtaCorrections \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC \
    --tests edu.duke.cs.osprey.packstar.TestPackStarConfig
sha256sum src/main/java/edu/duke/cs/osprey/packstar/PackStar*.java > "$RUN_ROOT/source.sha256"
echo "PASS: triple decomposition selector; reports in $RUN_ROOT/build"
