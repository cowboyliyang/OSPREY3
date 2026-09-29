#!/bin/bash
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -f "$RUN_ROOT/READY"
test -f "$BUILD_ROOT/MODEL_READY"
OUT="$RESULT_ROOT/shard_A$SLURM_ARRAY_JOB_ID/T$SLURM_ARRAY_TASK_ID"
mkdir -p "$OUT/tmp"
test ! -e "$OUT/STARTED"
touch "$OUT/STARTED"
export TMPDIR="$OUT/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "$BUILD_ROOT/edges.tsv" "$OUT/edges.tsv"
set +e
/usr/bin/time -o "$OUT/wall.time" -f 'elapsed=%e maxRssKiB=%M' \
 /home/users/lz280/java/jdk-17.0.2+8/bin/java -Xmx2g -XX:-UseSuperWord -XX:ActiveProcessorCount=4 \
 --add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED \
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
 "-Djava.io.tmpdir=$TMPDIR" "-Dcensus.cpus=$SLURM_CPUS_PER_TASK" \
 -cp "$(cat "$RUN_ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStarSparseNineCensus \
 shard "$BUILD_ROOT" "$OUT" "$BUILD_ROOT/config.tsv" "$SLURM_ARRAY_TASK_ID" "$BLOCK" > "$OUT/run.log" 2>&1
RC=$?
set -e
echo "$RC" > "$OUT/exit_code"
exit "$RC"
