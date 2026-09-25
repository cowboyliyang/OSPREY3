#!/bin/bash
# Reuse the already compiled dependencies; compile the final estimator and all
# focused test sources into a new immutable overlay before the comparison.
#SBATCH --job-name=packstar_forward_finalize
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${BASE_BUILD:?tested source build required}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_forward_20260913/$SLURM_JOB_ID
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/classes" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    "$REPO/src" "$REPO/lib" "$REPO/buildSrc" "$REPO/gradle" "$REPO/slurm" \
    "$REPO/gradlew" "$REPO/build.gradle.kts" "$REPO/settings.gradle.kts" "$REPO/docs" "$RUN_ROOT/source/"
CP=$(cat "$BASE_BUILD/test_classpath.txt")
# Recompile all Java sources into this snapshot. The base contributes compiled
# Kotlin and dependency jars; its standalone classpath-export invocation may
# have cleared Java outputs, so no Java class directory is assumed reusable.
cd "$RUN_ROOT/source"
rg --files src/main/java src/test/java -g '*.java' > "$RUN_ROOT/java-sources.txt"
"$JAVA_ROOT/bin/javac" -J-Xmx8g --add-modules jdk.incubator.foreign -cp "$CP" \
    -d "$RUN_ROOT/classes" @"$RUN_ROOT/java-sources.txt"
printf '%s:%s\n' "$RUN_ROOT/classes" "$CP" > "$RUN_ROOT/test_classpath.txt"
cd "$RUN_ROOT/source"
sha256sum src/main/java/edu/duke/cs/osprey/packstar/PackStar*.java > "$RUN_ROOT/source.sha256"
# Execute all 47 focused tests against the actual comparison classes.
"$JAVA_ROOT/bin/javac" -cp "$CP" -d "$RUN_ROOT/classes" \
    slurm/scripts/ForwardTripleTestRunner.java
"$JAVA_ROOT/bin/java" --add-modules jdk.incubator.foreign \
    -Djava.io.tmpdir="$TMPDIR" -cp "$(cat "$RUN_ROOT/test_classpath.txt")" ForwardTripleTestRunner
/bin/bash -n slurm/scripts/compare_forward_triples_20260913.slurm
/bin/bash -n slurm/scripts/summarize_forward_triples_20260913.slurm
python3 -m py_compile slurm/scripts/summarize_forward_triples_20260913.py
echo "PASS: final source overlay $RUN_ROOT"
