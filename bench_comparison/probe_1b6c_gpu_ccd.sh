#!/bin/bash
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3
OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CUDA_LIB=/usr/local/cuda-12.8/lib64
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx96g -Xms8g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
GEN=edu.duke.cs.osprey.markstar.bench.GenAuditConfs
MAIN=edu.duke.cs.osprey.markstar.bench.AuditBenchCpuGpu
mkdir -p "$LOGDIR"
./gradlew testClasses 2>&1 | tail -2
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.cp_1b6cprobe.txt"

# 1b6c (14 flex pos) -- the expensive case
B6_PDB=$PDBDIR/1b6c/1b6c.min.reduce.renum.pdb
B6_MUT='E965'
B6_FLEX='F994;F995;E903;E892;E957;E925;E963;E947;E922;E942;E940;E967;E948'
B6_CSV=$OUTDIR/audit_leaves/gen/1b6c_complex.csv
# 3u7y (5 flex pos) -- S11 anchor (~4.1x)
U7_PDB=$PDBDIR/3u7y/3u7y.min.reduce.renum.pdb
U7_MUT='G384;G382'
U7_FLEX='G385;L587;L649'
U7_CSV=$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv

WRAP="set -e; export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH; CP=\"\$(cat $LOGDIR/.cp_1b6cprobe.txt)\"; \
nvidia-smi --query-gpu=name,memory.total --format=csv || true; \
echo '##### GEN 1b6c confs #####'; mkdir -p $(dirname $B6_CSV); \
$JAVA $JARGS -Dosprey.bench.pdbPath=$B6_PDB -Dosprey.bench.designId=1b6c -Dosprey.bench.mutable='$B6_MUT' -Dosprey.bench.flexible='$B6_FLEX' -Dosprey.audit.state=complex -Dosprey.audit.output=$B6_CSV -Dosprey.gen.n=600 -cp \"\$CP\" $GEN; \
for DEV in cpu gpu; do \
  echo \"##### 1b6c \$DEV #####\"; \
  $JAVA $JARGS -Dosprey.bench.pdbPath=$B6_PDB -Dosprey.bench.designId=1b6c -Dosprey.bench.mutable='$B6_MUT' -Dosprey.bench.flexible='$B6_FLEX' -Dosprey.audit.state=complex -Dosprey.audit.input=$B6_CSV -Dosprey.bench.numCPUs=104 -Dosprey.audit.numGpus=8 -Dosprey.audit.streamsPerGpu=64 -Dosprey.audit.batchSize=1024 -Dosprey.audit.repeats=4 -Dosprey.audit.device=\$DEV -cp \"\$CP\" $MAIN || echo 1b6c_\${DEV}_FAILED; \
done; \
for DEV in cpu gpu; do \
  echo \"##### 3u7y \$DEV #####\"; \
  $JAVA $JARGS -Dosprey.bench.pdbPath=$U7_PDB -Dosprey.bench.designId=3u7y -Dosprey.bench.mutable='$U7_MUT' -Dosprey.bench.flexible='$U7_FLEX' -Dosprey.audit.state=Complex -Dosprey.audit.input=$U7_CSV -Dosprey.bench.numCPUs=104 -Dosprey.audit.numGpus=8 -Dosprey.audit.streamsPerGpu=64 -Dosprey.audit.batchSize=1024 -Dosprey.audit.maxConfs=600 -Dosprey.audit.repeats=4 -Dosprey.audit.device=\$DEV -cp \"\$CP\" $MAIN || echo 3u7y_\${DEV}_FAILED; \
done"

JID=$(sbatch --partition=grisman --account=grisman --nodes=1 --exclusive --mem=0 --constraint=a5000 \
  --time=02:00:00 --job-name=ccd_gpu_probe \
  --output=$LOGDIR/ccd_gpu_probe_%j.out --error=$LOGDIR/ccd_gpu_probe_%j.err \
  --wrap "$WRAP" | awk '{print $4}')
echo "Submitted CCD CPU/GPU probe: $JID -> $LOGDIR/ccd_gpu_probe_${JID}.out"
