#!/bin/bash
#SBATCH --job-name=ps_million_ccd
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -f "${BUILD_ROOT:?}/READY"
MODE=${MODE:-shard}
SYSTEM=${SYSTEM:-4wyu}
BLOCK=${BLOCK:-4096}
TASK=${SLURM_ARRAY_TASK_ID:-0}
OUT="${RESULT_ROOT:?}/${MODE}_A${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}/T$TASK"
mkdir -p "$OUT/tmp"
test ! -e "$OUT/STARTED"
touch "$OUT/STARTED"
export TMPDIR="$OUT/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "${BASH_SOURCE[0]}" "$OUT/runner.slurm"
scontrol show job "$SLURM_JOB_ID" > "$OUT/slurm.txt"
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ARGS=(--add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED
 -XX:-UseSuperWord -Xmx1g "-Djava.io.tmpdir=$TMPDIR" "-Dcensus.cpus=$SLURM_CPUS_PER_TASK"
 -cp "$(cat "$BUILD_ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStarMillionCensus
 "$MODE" "$BUILD_ROOT" "$OUT" "$BUILD_ROOT/screen/$SYSTEM.tsv" "$TASK" "$BLOCK")
if [[ -n "${RESUME_ROOT:-}" ]]; then ARGS+=("$RESUME_ROOT/T$TASK"); fi
printf '%q ' "$JAVA" "${ARGS[@]}" > "$OUT/command.sh"
set +e
/usr/bin/time -o "$OUT/wall.time" -f 'elapsed=%e maxRssKiB=%M' "$JAVA" "${ARGS[@]}" > "$OUT/run.log" 2>&1
RC=$?
set -e
printf '%s\n' "$RC" > "$OUT/exit_code"
exit "$RC"
