#!/bin/bash
#SBATCH --job-name=ps-sample-accounting
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_sample_accounting_20261003/check_${SLURM_JOB_ID}
# Source/fixture copy ~130 MiB / several thousand files; build/reports <2 GiB expected.
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
cd "$REPO"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib examples buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
git diff --binary HEAD > "$RUN_ROOT/source.patch"
git ls-files --others --exclude-standard > "$RUN_ROOT/untracked_sources.txt"
cd "$RUN_ROOT/source"
./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    -I slurm/scripts/gradle_xtmp.init.gradle -I slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" -DtestMaxHeap=12g \
    --tests 'edu.duke.cs.osprey.packstar.TestPackStar*'
echo "SAMPLE_ACCOUNTING_VALIDATION=$RUN_ROOT"
