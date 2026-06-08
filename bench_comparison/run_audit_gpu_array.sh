#!/bin/bash
# Run S11 CCD GPU audit as a Slurm array: one shard per GPU, then merge outputs.
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CUDA_LIB=/usr/local/cuda-12.8/lib64
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx16g -Xms1g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
AUDIT_MAIN=edu.duke.cs.osprey.markstar.bench.AuditLeafCCD

NUM_SHARDS=${NUM_SHARDS:-48}
MAX_CONCURRENT=${MAX_CONCURRENT:-$NUM_SHARDS}
ARRAY_SPEC=${ARRAY_SPEC:-0-$((NUM_SHARDS - 1))%$MAX_CONCURRENT}
STREAMS=${STREAMS:-64}
BATCH_SIZE=${BATCH_SIZE:-1024}
WARMUP_CONFS=${WARMUP_CONFS:-8}
CPUS=${CPUS:-8}
MEM=${MEM:-20G}
TIME_LIMIT=${TIME_LIMIT:-00:30:00}
NODELIST=${NODELIST:-fennario-[01-06]}
RUN_TAG=${RUN_TAG:-3u7y_Complex_seq00034_fennario${NUM_SHARDS}_$(date +%Y%m%d_%H%M%S)}

DID=3u7y
STATE=Complex
PDB=$PDBDIR/$DID/${DID}.min.reduce.renum.pdb
INPUT=$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv
MUT='G384;G382'
FLEX='G385;L587;L649'
EXPECTED_ROWS=${EXPECTED_ROWS:-2752}

mkdir -p "$LOGDIR" "$OUTDIR/audit_results"
SHARDDIR=$OUTDIR/audit_results/${RUN_TAG}_shards
MERGED=$OUTDIR/audit_results/${RUN_TAG}_merged.csv
mkdir -p "$SHARDDIR"

./gradlew testClasses 2>&1 | tail -2

CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
CPFILE=$LOGDIR/.classpath_audit_array_${RUN_TAG}.txt
echo "$CP" > "$CPFILE"

ARRAY_WRAP="set -eu; \
export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH; \
TASK_TMP=$OUTDIR/tmp/audit_gpu_array_${RUN_TAG}_\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}; \
mkdir -p \"\$TASK_TMP/cuda_cache\"; \
export TMPDIR=\"\$TASK_TMP\"; \
export CUDA_CACHE_PATH=\"\$TASK_TMP/cuda_cache\"; \
CP=\"\$(cat $CPFILE)\"; \
OUT=$SHARDDIR/shard_\${SLURM_ARRAY_TASK_ID}.csv; \
echo \"##### ARRAY SHARD \${SLURM_ARRAY_TASK_ID}/$NUM_SHARDS host=\$(hostname) streams=$STREAMS batchSize=$BATCH_SIZE warmupConfs=$WARMUP_CONFS tmp=\$TASK_TMP output=\$OUT #####\"; \
nvidia-smi --query-gpu=name,index,driver_version,memory.total --format=csv || true; \
/usr/bin/time -p $JAVA $JARGS \"-Djava.io.tmpdir=\$TASK_TMP\" \"-XX:ErrorFile=\$TASK_TMP/hs_err_pid%p.log\" \
\"-Dosprey.bench.pdbPath=$PDB\" \
\"-Dosprey.bench.mutable=$MUT\" \
\"-Dosprey.bench.flexible=$FLEX\" \
\"-Dosprey.bench.designId=$DID\" \
\"-Dosprey.audit.state=$STATE\" \
\"-Dosprey.audit.input=$INPUT\" \
\"-Dosprey.bench.numCPUs=$CPUS\" \
\"-Dosprey.audit.output=\$OUT\" \
\"-Dosprey.audit.device=gpu\" \
\"-Dosprey.audit.streamsPerGpu=$STREAMS\" \
\"-Dosprey.audit.batchSize=$BATCH_SIZE\" \
\"-Dosprey.audit.warmupConfs=$WARMUP_CONFS\" \
\"-Dosprey.audit.numShards=$NUM_SHARDS\" \
\"-Dosprey.audit.shardIndex=\${SLURM_ARRAY_TASK_ID}\" \
-cp \"\$CP\" $AUDIT_MAIN; \
wc -l \"\$OUT\""

ARRAY_JID=$(sbatch --parsable --partition=grisman --account=grisman \
    --array="$ARRAY_SPEC" \
    --nodes=1 --nodelist="$NODELIST" --gres=gpu:a5000:1 \
    --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME_LIMIT" \
    --job-name=audit_gpu_arr \
    --output=$LOGDIR/audit_gpu_array_${RUN_TAG}_%A_%a.out \
    --error=$LOGDIR/audit_gpu_array_${RUN_TAG}_%A_%a.err \
    --wrap "$ARRAY_WRAP")

MERGE_WRAP="set -eu; \
MERGED=$MERGED; \
SHARDDIR=$SHARDDIR; \
NUM_SHARDS=$NUM_SHARDS; \
EXPECTED_ROWS=$EXPECTED_ROWS; \
echo \"##### MERGE \$SHARDDIR -> \$MERGED #####\"; \
found=\$(find \"\$SHARDDIR\" -maxdepth 1 -name 'shard_*.csv' | wc -l); \
echo \"found_shards=\$found expected_shards=\$NUM_SHARDS\"; \
test \"\$found\" -eq \"\$NUM_SHARDS\"; \
first=1; \
rm -f \"\$MERGED\"; \
for f in \$(find \"\$SHARDDIR\" -maxdepth 1 -name 'shard_*.csv' | sort -V); do \
  if [ \$first -eq 1 ]; then cat \"\$f\" > \"\$MERGED\"; first=0; else tail -n +2 \"\$f\" >> \"\$MERGED\"; fi; \
done; \
wc -l \"\$MERGED\"; \
data_rows=\$(( \$(wc -l < \"\$MERGED\") - 1 )); \
echo \"data_rows=\$data_rows expected_rows=\$EXPECTED_ROWS\"; \
test \"\$data_rows\" -eq \"\$EXPECTED_ROWS\"; \
awk -F, 'NR>1 {count[\$15]++} END {for (s in count) print s,count[s]}' \"\$MERGED\""

MERGE_JID=$(sbatch --parsable --partition=grisman --account=grisman \
    --dependency=afterok:$ARRAY_JID --nodes=1 --cpus-per-task=2 --mem=4G --time=00:10:00 \
    --job-name=audit_gpu_merge \
    --output=$LOGDIR/audit_gpu_array_${RUN_TAG}_merge_%j.out \
    --error=$LOGDIR/audit_gpu_array_${RUN_TAG}_merge_%j.err \
    --wrap "$MERGE_WRAP")

echo "Submitted GPU audit array: $ARRAY_JID"
echo "  run_tag: $RUN_TAG"
echo "  array:   $ARRAY_SPEC"
echo "  shards:  $SHARDDIR"
echo "  logs:    $LOGDIR/audit_gpu_array_${RUN_TAG}_${ARRAY_JID}_*.out"
echo "Submitted merge job: $MERGE_JID"
echo "  merged:  $MERGED"
