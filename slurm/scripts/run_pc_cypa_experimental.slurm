#!/bin/bash
#SBATCH --job-name=pc-cypa-experimental
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
PREP=/usr/xtmp/lz280/packstar_functional_validation/experimental_12552449
RUN_ROOT=/usr/xtmp/lz280/packstar_functional_validation/cypa_${SLURM_JOB_ID}
SOURCE=$PREP/pdbs_prepped/3k0n_m001
test -s "$SOURCE/3k0n_m001.min.reduce.renum.pdb"
mkdir -p "$RUN_ROOT"/{input,output,source_snapshot,gradle_home/caches,tmp,frequency_severity}
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="$RUN_ROOT/gradle_home"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
export PACKSTAR_FORMAL_CLASSPATH_FILE="$RUN_ROOT/classpath.txt"
cp -a /usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache/wrapper "$GRADLE_USER_HOME/"
cp -a /usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache/caches/modules-2 "$GRADLE_USER_HOME/caches/"
cp "$SOURCE/3k0n_m001.min.reduce.renum.pdb" "$RUN_ROOT/input/cypa.pdb"
cp "$SOURCE/3k0n_m001_strip_reduce_prep_renum.txt" "$RUN_ROOT/input/renum.txt"
cp "$PREP/sources.tsv" "$RUN_ROOT/input/sources.tsv"
cp "$REPO/docs/packstar_pc_cypa_experimental_v1.json" "$RUN_ROOT/protocol.json"
cd "$REPO"
# Entire source is rebuilt because Java may inline temperature constants.
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ --exclude=/.gradle/ --exclude='*.class' \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source_snapshot/"
git diff --binary > "$RUN_ROOT/tracked_changes.patch"
git status --short > "$RUN_ROOT/source_status.txt"
git rev-parse HEAD > "$RUN_ROOT/git_head.txt"
cd "$RUN_ROOT/source_snapshot"
mkdir -p examples/python.KStar
cp "$REPO/examples/python.KStar/2RL0.min.reduce.pdb" examples/python.KStar/
THERMAL=src/main/java/edu/duke/cs/osprey/kstar/pfunc/BoltzmannCalculator.java
cp "$THERMAL" "$RUN_ROOT/BoltzmannCalculator.original.java"
sed -i 's/public static final double TClassic = TRoom;/public static final double TClassic = 283.15;/' "$THERMAL"
for THERMAL_AUX in src/main/java/edu/duke/cs/osprey/branchdp/BranchDPScorer.java src/main/java/edu/duke/cs/osprey/branchdp/DecompSearchNode.java; do
    cp "$THERMAL_AUX" "$RUN_ROOT/$(basename "$THERMAL_AUX").original"
    sed -i 's/BoltzmannCalculator.RClassic \* 298.0/BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic/g' "$THERMAL_AUX"
done
diff -u "$RUN_ROOT/BoltzmannCalculator.original.java" "$THERMAL" > "$RUN_ROOT/temperature_overlay.patch" || test "$?" = 1
find src slurm -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_ROOT/source_sha256.txt"
sha256sum "$RUN_ROOT/input/cypa.pdb" "$RUN_ROOT/input/renum.txt" "$RUN_ROOT/protocol.json" > "$RUN_ROOT/input_sha256.txt"
./gradlew test --offline --no-daemon --max-workers=1 \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalGeometry \
    --tests edu.duke.cs.osprey.packstar.TestPackStarResult.kstarRetainsFunctionalEstimateAndBackendCounters \
    writeFormalCaseClasspath --project-cache-dir "$RUN_ROOT/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/packstar_formal_classpath.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -DtestMaxHeap=8g -Dorg.gradle.jvmargs="-Xmx12g -Djava.io.tmpdir=$TMPDIR" \
    > "$RUN_ROOT/build.stdout" 2> "$RUN_ROOT/build.stderr"
CP=$(cat "$RUN_ROOT/classpath.txt")
/usr/bin/time -v "$JAVA_HOME/bin/java" \
    --add-opens java.base/java.util=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
    -Xmx48g -Djava.io.tmpdir="$TMPDIR" -Dcypa.threads=8 \
    -Dpackstar.pac.samples=2048 -Dpackstar.pac.maxEstSamples=4000 \
    -Dpackstar.pac.unreachableCap=400 -Dpackstar.pac.confidence=0.05 \
    -Dpackstar.pac.targetEpsilon=0.2 -Dpackstar.pac.randomSeed=42 \
    -Dpackstar.pac.sampling.gpu=false -Dpackstar.pac.sampling.progress=false \
    -Dpackstar.pac.sampling.threads=8 -Dbranchdp.dp.gpu=false \
    -Dpackstar.pac.frequencySeverity.outputDir="$RUN_ROOT/frequency_severity" \
    -cp "$CP" edu.duke.cs.osprey.packstar.PackStarCypAExperimentalCase \
    "$RUN_ROOT/input/cypa.pdb" "$RUN_ROOT/input/renum.txt" "$RUN_ROOT/output" \
    > "$RUN_ROOT/runner.stdout" 2> "$RUN_ROOT/runner.stderr"
test -s "$RUN_ROOT/output/COMPUTATION_FINISHED"
