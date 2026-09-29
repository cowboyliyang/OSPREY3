#!/bin/bash
# Source snapshot: approximately 100 MiB / several thousand files; build/cache on xtmp.
#SBATCH --job-name=ps-auto-check
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_auto_preflight/check_${SLURM_JOB_ID}
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
export PACKSTAR_CLASSPATH_FILE="$RUN_ROOT/test_classpath.txt"
cd "$REPO"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
cd "$RUN_ROOT/source"
bash -n slurm/scripts/run_packstar_auto_preflight.slurm
./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    -I slurm/scripts/gradle_xtmp.init.gradle -I slurm/scripts/verify_ccd.init.gradle \
    -I slurm/h200/classpath.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" -DtestMaxHeap=8g \
    --tests edu.duke.cs.osprey.packstar.TestPackStarAutoPreflight \
    --tests edu.duke.cs.osprey.packstar.TestPackStarAdmissionDecision \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDpAdmission \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleDecompositionCosts \
    writeExperimentClasspath
echo "AUTO_PREFLIGHT_CLASSPATH=$PACKSTAR_CLASSPATH_FILE"
ln -sfn "$RUN_ROOT" /usr/xtmp/lz280/packstar_auto_preflight/current
