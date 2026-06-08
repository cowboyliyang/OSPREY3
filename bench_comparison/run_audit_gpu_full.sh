#!/bin/bash
# Run a tuned S11 CCD GPU audit on one audit-leaf CSV.
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CUDA_LIB=/usr/local/cuda-12.8/lib64
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx32g -Xms2g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
BENCH_MAIN=edu.duke.cs.osprey.markstar.bench.AuditBenchCpuGpu
AUDIT_MAIN=edu.duke.cs.osprey.markstar.bench.AuditLeafCCD

FULL_STREAMS=${FULL_STREAMS:-64}
FULL_BATCH_SIZE=${FULL_BATCH_SIZE:-1024}
WARMUP_CONFS=${WARMUP_CONFS:-64}
SMOKE_CONFS=${SMOKE_CONFS:-8}
SMOKE_REPEATS=${SMOKE_REPEATS:-1}

DID=3u7y
STATE=Complex
PDB=$PDBDIR/$DID/${DID}.min.reduce.renum.pdb
INPUT=$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv
MUT='G384;G382'
FLEX='G385;L587;L649'

mkdir -p "$LOGDIR" "$OUTDIR/audit_results"
./gradlew testClasses 2>&1 | tail -2

CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_auditbench.txt"

COMMON="-Dosprey.bench.pdbPath=$PDB -Dosprey.bench.mutable='$MUT' -Dosprey.bench.flexible='$FLEX' -Dosprey.bench.designId=$DID -Dosprey.audit.state=$STATE -Dosprey.audit.input=$INPUT -Dosprey.bench.numCPUs=8"

WRAP="set -eu; \
export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH; \
CP=\"\$(cat $LOGDIR/.classpath_auditbench.txt)\"; \
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv || true; \
echo \"##### GPU SMOKE streams=$FULL_STREAMS batchSize=$FULL_BATCH_SIZE maxConfs=$SMOKE_CONFS #####\"; \
$JAVA $JARGS $COMMON -Dosprey.audit.device=gpu -Dosprey.audit.streamsPerGpu=$FULL_STREAMS -Dosprey.audit.batchSize=$FULL_BATCH_SIZE -Dosprey.audit.maxConfs=$SMOKE_CONFS -Dosprey.audit.repeats=$SMOKE_REPEATS -cp \"\$CP\" $BENCH_MAIN; \
FULL_OUT=$OUTDIR/audit_results/gpu_full_${DID}_${STATE}_seq00034_\${SLURM_JOB_ID}.csv; \
echo \"##### FULL AUDITLEAFCCD GPU streams=$FULL_STREAMS batchSize=$FULL_BATCH_SIZE warmupConfs=$WARMUP_CONFS output=\$FULL_OUT #####\"; \
/usr/bin/time -p $JAVA $JARGS $COMMON -Dosprey.audit.output=\$FULL_OUT -Dosprey.audit.device=gpu -Dosprey.audit.streamsPerGpu=$FULL_STREAMS -Dosprey.audit.batchSize=$FULL_BATCH_SIZE -Dosprey.audit.warmupConfs=$WARMUP_CONFS -cp \"\$CP\" $AUDIT_MAIN; \
echo \"##### FULL OUTPUT ROWS #####\"; \
wc -l \$FULL_OUT"

JID=$(sbatch --partition=grisman --account=grisman --nodes=1 --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=01:00:00 \
    --nodelist=fennario-01,fennario-02,fennario-03 --job-name=audit_gpu_full \
    --output=$LOGDIR/audit_gpu_full_%j.out --error=$LOGDIR/audit_gpu_full_%j.err \
    --wrap "$WRAP" | awk '{print $4}')
echo "Submitted full GPU audit: $JID -> $LOGDIR/audit_gpu_full_${JID}.out"
