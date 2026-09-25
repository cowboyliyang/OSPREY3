#!/bin/bash
#SBATCH --job-name=pc-cypa-resume
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
BUNDLE=/usr/xtmp/lz280/packstar_functional_validation/cypa_12552499
RUN_ROOT=/usr/xtmp/lz280/packstar_functional_validation/cypa_${SLURM_JOB_ID}
mkdir -p "$RUN_ROOT"/{tmp,output,frequency_severity}
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="$BUNDLE/gradle_home"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUNDLE/build"
export PACKSTAR_FORMAL_CLASSPATH_FILE="$RUN_ROOT/classpath.txt"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/resume_pc_cypa_experimental.slurm "$RUN_ROOT/submission.slurm"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_pc_cypa.py "$RUN_ROOT/summarize.py"
cp "$BUNDLE/protocol.json" "$RUN_ROOT/protocol.json"
sed -i 's/"epsilon": 0.1/"epsilon": 0.2/' "$RUN_ROOT/protocol.json"
printf 'source_bundle\t%s\nreason\tExclude unrelated 4WYU 298.15K-specific protocol fixture from 283.15K experiment\n' "$BUNDLE" > "$RUN_ROOT/execution.tsv"
cd "$BUNDLE/source_snapshot"
sha256sum -c "$BUNDLE/source_sha256.txt" > "$RUN_ROOT/source_verification.txt"
sha256sum -c "$BUNDLE/input_sha256.txt" > "$RUN_ROOT/input_verification.txt"
./gradlew test --offline --no-daemon --max-workers=1 \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalGeometry \
    --tests edu.duke.cs.osprey.packstar.TestPackStarResult.kstarRetainsFunctionalEstimateAndBackendCounters \
    writeFormalCaseClasspath --project-cache-dir "$BUNDLE/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/packstar_formal_classpath.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -DtestMaxHeap=8g -Dorg.gradle.jvmargs="-Xmx12g -Djava.io.tmpdir=$TMPDIR" \
    > "$RUN_ROOT/build.stdout" 2> "$RUN_ROOT/build.stderr"
cp -a "$BUNDLE/build/osprey3-root/test-results" "$RUN_ROOT/regression-results"
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
    "$BUNDLE/input/cypa.pdb" "$BUNDLE/input/renum.txt" "$RUN_ROOT/output" \
    > "$RUN_ROOT/runner.stdout" 2> "$RUN_ROOT/runner.stderr"
/home/users/lz280/miniconda3/envs/confdiff/bin/python "$RUN_ROOT/summarize.py" "$RUN_ROOT" > "$RUN_ROOT/summary.stdout"
