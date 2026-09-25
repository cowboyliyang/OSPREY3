#!/bin/bash
#SBATCH --job-name=ccd3cal_trp
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -f "${BUILD_ROOT:?}/READY"
TASK=${SLURM_ARRAY_TASK_ID:?}
OUT="${RESULT_ROOT:?}/shard_A${SLURM_ARRAY_JOB_ID}/T$TASK"
mkdir -p "$OUT/tmp"
test ! -e "$OUT/STARTED"
touch "$OUT/STARTED"
export TMPDIR="$OUT/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "${BASH_SOURCE[0]}" "$OUT/runner.slurm"
scontrol show job "$SLURM_JOB_ID" > "$OUT/slurm.txt"
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ARGS=(--add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED
 -XX:-UseSuperWord -Xmx2g -XX:ActiveProcessorCount=4 "-Djava.io.tmpdir=$TMPDIR" "-Dcensus.cpus=$SLURM_CPUS_PER_TASK" "-Dcensus.target=${TARGET:?}"
 -cp "$(cat "$BUILD_ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStar3calProteinCensus
 shard "$BUILD_ROOT" "$OUT" "$TASK" "${BLOCK:-32768}")
if [[ -n "${RESUME_ROOT:-}" ]]; then ARGS+=("$RESUME_ROOT/T$TASK"); fi
printf '%q ' "$JAVA" "${ARGS[@]}" > "$OUT/command.sh"
set +e
/usr/bin/time -o "$OUT/wall.time" -f 'elapsed=%e maxRssKiB=%M' "$JAVA" "${ARGS[@]}" > "$OUT/run.log" 2>&1
RC=$?
set -e
printf '%s\n' "$RC" > "$OUT/exit_code"
exit "$RC"
