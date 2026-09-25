#!/bin/bash
#SBATCH --job-name=pc-cypa-array
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail

test -n "${SLURM_JOB_ID:-}"
: "${BUILD_ROOT:?BUILD_ROOT must be exported}"
test -s "$BUILD_ROOT/classpath.txt"
test -s "$BUILD_ROOT/build_manifest.tsv"

ARRAY_ROOT=/usr/xtmp/lz280/packstar_functional_validation/cypa_parallel_${SLURM_ARRAY_JOB_ID}
RUN_ROOT="$ARRAY_ROOT/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$RUN_ROOT"/{input,output,frequency_severity,tmp}
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

cp "$BUILD_ROOT/protocol.json" "$RUN_ROOT/protocol.json"
cp "$BUILD_ROOT/input/cypa.pdb" "$RUN_ROOT/input/cypa.pdb"
cp "$BUILD_ROOT/input/renum.txt" "$RUN_ROOT/input/renum.txt"
cp "$BUILD_ROOT/input/sources.tsv" "$RUN_ROOT/input/sources.tsv"
cp "$BUILD_ROOT/classpath.txt" "$RUN_ROOT/classpath.txt"
printf 'array_job\t%s\narray_task\t%s\nvariant_index\t%s\nbuild_root\t%s\ncpus\t%s\ntarget_epsilon\t0.2\n' \
    "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$SLURM_ARRAY_TASK_ID" "$BUILD_ROOT" \
    "${SLURM_CPUS_PER_TASK:-48}" > "$RUN_ROOT/task_manifest.tsv"

CP=$(cat "$RUN_ROOT/classpath.txt")
/usr/bin/time -v "$JAVA_HOME/bin/java" \
    --add-opens java.base/java.util=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
    -Xmx48g -Djava.io.tmpdir="$TMPDIR" \
    -Dcypa.threads="${SLURM_CPUS_PER_TASK:-48}" \
    -Dpackstar.pac.samples=2048 -Dpackstar.pac.maxEstSamples=4000 \
    -Dpackstar.pac.unreachableCap=400 -Dpackstar.pac.confidence=0.05 \
    -Dpackstar.pac.targetEpsilon=0.2 -Dpackstar.pac.randomSeed=42 \
    -Dpackstar.pac.sampling.gpu=false -Dpackstar.pac.sampling.progress=false \
    -Dpackstar.pac.sampling.threads="${SLURM_CPUS_PER_TASK:-48}" \
    -Dbranchdp.dp.gpu=false \
    -Dpackstar.pac.frequencySeverity.outputDir="$RUN_ROOT/frequency_severity" \
    -cp "$CP" edu.duke.cs.osprey.packstar.PackStarCypAExperimentalCase \
    "$RUN_ROOT/input/cypa.pdb" "$RUN_ROOT/input/renum.txt" "$RUN_ROOT/output" \
    "$SLURM_ARRAY_TASK_ID" \
    > "$RUN_ROOT/runner.stdout" 2> "$RUN_ROOT/runner.stderr"
test -s "$RUN_ROOT/output/COMPUTATION_FINISHED"
