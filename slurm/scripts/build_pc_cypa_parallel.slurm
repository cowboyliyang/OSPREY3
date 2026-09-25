#!/bin/bash
#SBATCH --job-name=pc-cypa-build
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail

test -n "${SLURM_JOB_ID:-}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
PREP=/usr/xtmp/lz280/packstar_functional_validation/experimental_12552449
SOURCE="$PREP/pdbs_prepped/3k0n_m001"
ROOT=/usr/xtmp/lz280/packstar_functional_validation/cypa_parallel_build_${SLURM_JOB_ID}
CACHE=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache

test -s "$SOURCE/3k0n_m001.min.reduce.renum.pdb"
test -s "$SOURCE/3k0n_m001_strip_reduce_prep_renum.txt"
mkdir -p "$ROOT"/{input,source_snapshot,gradle_home/caches,tmp,frequency_severity}
export TMPDIR="$ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="$ROOT/gradle_home"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$ROOT/build"
export PACKSTAR_FORMAL_CLASSPATH_FILE="$ROOT/classpath.txt"

cp -a "$CACHE/wrapper" "$GRADLE_USER_HOME/"
cp -a "$CACHE/caches/modules-2" "$GRADLE_USER_HOME/caches/"
cp "$SOURCE/3k0n_m001.min.reduce.renum.pdb" "$ROOT/input/cypa.pdb"
cp "$SOURCE/3k0n_m001_strip_reduce_prep_renum.txt" "$ROOT/input/renum.txt"
cp "$PREP/sources.tsv" "$ROOT/input/sources.tsv"
cp "$REPO/docs/packstar_pc_cypa_experimental_v1.json" "$ROOT/protocol.json"
cp "$REPO/slurm/scripts/build_pc_cypa_parallel.slurm" "$ROOT/submission.slurm"

cd "$REPO"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ --exclude=/.gradle/ --exclude='*.class' \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$ROOT/source_snapshot/"
git diff --binary > "$ROOT/tracked_changes.patch"
git status --short > "$ROOT/source_status.txt"
git rev-parse HEAD > "$ROOT/git_head.txt"
cd "$ROOT/source_snapshot"
mkdir -p examples/python.KStar
cp "$REPO/examples/python.KStar/2RL0.min.reduce.pdb" examples/python.KStar/

# The CypA experiment uses 283.15 K.  These overlays make the temperature
# consistent in both Boltzmann factors and branch-and-bound scoring.
THERMAL=src/main/java/edu/duke/cs/osprey/kstar/pfunc/BoltzmannCalculator.java
cp "$THERMAL" "$ROOT/BoltzmannCalculator.original.java"
sed -i 's/public static final double TClassic = TRoom;/public static final double TClassic = 283.15;/' "$THERMAL"
for THERMAL_AUX in \
    src/main/java/edu/duke/cs/osprey/branchdp/BranchDPScorer.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DecompSearchNode.java; do
    cp "$THERMAL_AUX" "$ROOT/$(basename "$THERMAL_AUX").original"
    sed -i 's/BoltzmannCalculator.RClassic \* 298.0/BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic/g' "$THERMAL_AUX"
done
diff -u "$ROOT/BoltzmannCalculator.original.java" "$THERMAL" > "$ROOT/temperature_overlay.patch" || test "$?" = 1
find src slurm -type f -print0 | sort -z | xargs -0 sha256sum > "$ROOT/source_sha256.txt"
sha256sum "$ROOT/input/cypa.pdb" "$ROOT/input/renum.txt" "$ROOT/protocol.json" > "$ROOT/input_sha256.txt"

./gradlew test --offline --no-daemon --max-workers=1 \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalGeometry \
    --tests edu.duke.cs.osprey.packstar.TestPackStarResult.kstarRetainsFunctionalEstimateAndBackendCounters \
    writeFormalCaseClasspath --project-cache-dir "$ROOT/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/packstar_formal_classpath.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -DtestMaxHeap=8g -Dorg.gradle.jvmargs="-Xmx12g -Djava.io.tmpdir=$TMPDIR" \
    > "$ROOT/build.stdout" 2> "$ROOT/build.stderr"
test -s "$ROOT/classpath.txt"
printf 'build_root\t%s\nsource_variant_selector\ttrue\ntarget_epsilon\t0.2\n' "$ROOT" > "$ROOT/build_manifest.tsv"
