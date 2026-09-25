#!/bin/bash
# Called only inside Slurm; reuse the validated source/build/input bundle.
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
BUNDLE=$1
RUN_ROOT=$2
JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
mkdir -p "$RUN_ROOT/output" "$RUN_ROOT/tmp" "$RUN_ROOT/frequency_severity"
export TMPDIR="$RUN_ROOT/tmp"
cd "$BUNDLE/source_snapshot"
sha256sum -c "$BUNDLE/source_sha256.txt" > "$RUN_ROOT/source_verification.txt"
CP=$(cat "$BUNDLE/classpath.txt")
PREFLIGHT=false
if [ "${VALIDATE_ONLY:-0}" = 1 ]; then PREFLIGHT=true; fi
printf 'job_id\t%s\nbundle\t%s\npreflight\t%s\n' "$SLURM_JOB_ID" "$BUNDLE" "$PREFLIGHT" > "$RUN_ROOT/execution.tsv"
/usr/bin/time -v "$JAVA_HOME/bin/java" \
    --add-opens java.base/java.util=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
    -Xmx112g -Xms8g -Djava.io.tmpdir="$TMPDIR" \
    -Dpackstar.pac.samples=1000 -Dpackstar.pac.confidence=0.05 \
    -Dpackstar.pac.targetEpsilon=0.1 -Dpackstar.pac.randomSeed=42 \
    -Dpackstar.pac.sampling.gpu=false -Dpackstar.pac.sampling.progress=false \
    -Dpackstar.pac.sampling.threads="$SLURM_CPUS_PER_TASK" \
    -Dpackstar.pac.frequencySeverity.outputDir="$RUN_ROOT/frequency_severity" \
    -Dpc4wyu.proteinPdb="$BUNDLE/input/protein_A.pdb" \
    -Dpc4wyu.ligandPdb="$BUNDLE/input/peptide_D.pdb" \
    -Dpc4wyu.protocol="$BUNDLE/protocol.json" \
    -Dpc4wyu.preflight="$PREFLIGHT" -Dpc4wyu.output="$RUN_ROOT/output" \
    -Dpc4wyu.threads="$SLURM_CPUS_PER_TASK" -Dpc4wyu.epsilon=0.1 \
    -Dpc4wyu.maxNumConfs=2147483647 -cp "$CP" \
    edu.duke.cs.osprey.packstar.PackStar4wyuFormalCase \
    > "$RUN_ROOT/runner.stdout" 2> "$RUN_ROOT/runner.stderr"
if [ "$PREFLIGHT" = true ]; then
    printf 'Regression tests and 4WYU input/space preflight passed\n' > "$RUN_ROOT/VALIDATED"
else
    test -f "$RUN_ROOT/output/COMPLETE"
    test -s "$RUN_ROOT/output/validation.tsv"
    test -s "$RUN_ROOT/output/selections.tsv"
    echo "Execution complete; inspect validation.tsv before drawing scientific conclusions."
fi
