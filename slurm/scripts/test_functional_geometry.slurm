#!/bin/bash
# Pre-freeze geometry/reference regression; no formal design benchmark.
#SBATCH --job-name=packstar_pc_geometry
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
RUN_ROOT=/usr/xtmp/lz280/packstar_functional_validation/geometry_$SLURM_JOB_ID
mkdir -p "$RUN_ROOT/source/examples/python.KStar" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
cd "$REPO"
# Same bounded source snapshot as the selector regression (~60 MiB), plus
# one existing <1 MiB PDB. Allow 2 GiB / 20000 files for source and build outputs.
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
cp examples/python.KStar/2RL0.min.reduce.pdb "$RUN_ROOT/source/examples/python.KStar/"
git diff --binary > "$RUN_ROOT/tracked_changes.patch"
git status --short > "$RUN_ROOT/source_status.txt"
cd "$RUN_ROOT/source"
find src -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_ROOT/source_sha256.txt"
sha256sum examples/python.KStar/2RL0.min.reduce.pdb > "$RUN_ROOT/input_sha256.txt"
date --iso-8601=seconds
hostname
/usr/bin/time -v ./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    -DtestMaxHeap=8g \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalGeometry \
    --tests edu.duke.cs.osprey.packstar.TestPackStarConformationDedup
echo "PASS: P_C geometry/reference regression; reports in $RUN_ROOT/build"
