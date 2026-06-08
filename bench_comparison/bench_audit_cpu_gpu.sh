#!/bin/bash
# CPU-vs-GPU CCD throughput probe on an existing S11 audit-leaf CSV.
# Runs (same node, identical hardware): GPU smoke -> CPU full -> GPU full @ useful stream counts.
# Each config is its own JVM with '|| true' so a GPU failure cannot abort the CPU run.
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CUDA_LIB=/usr/local/cuda-12.8/lib64
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx32g -Xms2g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.bench.AuditBenchCpuGpu
CPU_REPEATS=${CPU_REPEATS:-2}
GPU16_REPEATS=${GPU16_REPEATS:-2}
GPU64_REPEATS=${GPU64_REPEATS:-2}

# 3u7y / Complex, largest audit-leaf CSV (~2752 confs)
DID=3u7y
PDB=$PDBDIR/$DID/${DID}.min.reduce.renum.pdb
INPUT=$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv
MUT='G384;G382'
FLEX='G385;L587;L649'

mkdir -p "$LOGDIR"
./gradlew testClasses 2>&1 | tail -2
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_auditbench.txt"

COMMON="-Dosprey.bench.pdbPath=$PDB -Dosprey.bench.mutable='$MUT' -Dosprey.bench.flexible='$FLEX' -Dosprey.bench.designId=$DID -Dosprey.audit.state=Complex -Dosprey.audit.input=$INPUT -Dosprey.bench.numCPUs=8"

WRAP="export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && CP=\"\$(cat $LOGDIR/.classpath_auditbench.txt)\" && nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv || true && \
echo '##### [1] GPU SMOKE (8 confs) #####' && $JAVA $JARGS $COMMON -Dosprey.audit.device=gpu -Dosprey.audit.streamsPerGpu=16 -Dosprey.audit.maxConfs=8 -Dosprey.audit.repeats=1 -cp \"\$CP\" $MAIN || echo GPU_SMOKE_FAILED && \
echo '##### [2] CPU full (8 threads) #####' && $JAVA $JARGS $COMMON -Dosprey.audit.device=cpu -Dosprey.audit.repeats=$CPU_REPEATS -cp \"\$CP\" $MAIN || echo CPU_FAILED && \
echo '##### [3] GPU full streams=16 #####' && $JAVA $JARGS $COMMON -Dosprey.audit.device=gpu -Dosprey.audit.streamsPerGpu=16 -Dosprey.audit.repeats=$GPU16_REPEATS -cp \"\$CP\" $MAIN || echo GPU16_FAILED && \
echo '##### [4] GPU full streams=64 batchSize=1024 #####' && $JAVA $JARGS $COMMON -Dosprey.audit.device=gpu -Dosprey.audit.streamsPerGpu=64 -Dosprey.audit.batchSize=1024 -Dosprey.audit.repeats=$GPU64_REPEATS -cp \"\$CP\" $MAIN || echo GPU64_FAILED"

JID=$(sbatch --partition=grisman --account=grisman --nodes=1 --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=02:00:00 \
    --nodelist=fennario-01,fennario-02,fennario-03 --job-name=audit_cpugpu \
    --output=$LOGDIR/audit_cpugpu_%j.out --error=$LOGDIR/audit_cpugpu_%j.err \
    --wrap "$WRAP" | awk '{print $4}')
echo "Submitted audit CPU/GPU probe: $JID -> $LOGDIR/audit_cpugpu_${JID}.out"
