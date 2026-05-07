#!/bin/bash
# Submit GNN benchmark: CCD baseline vs Strategy7 (decoupled GPU batch)
# Uses A5000 on grisman partition.
# Usage: bash bench_s7.sh [gpuBatchSize1 gpuBatchSize2 ...]
# Default batch sizes: 500, 1000, 2000

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

SLURM="--partition=grisman --account=grisman --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=02:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

PMODEL="gnn_data/2RL0_all20_4pos_merged/protein/model/gnn_model.onnx"
CMODEL="gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_model.onnx"
NUMSEQS=5
CONFSPACE=highrot
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx32g -Xms2g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.RunBenchmark

# GPU batch sizes from args, or defaults
if [ $# -gt 0 ]; then
    BATCH_SIZES=("$@")
else
    BATCH_SIZES=(500 1000 2000)
fi

# Build
echo "Building..."
./gradlew testClasses 2>&1 | tail -3

# Classpath
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do
    CP="$CP:$jar"
done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do
    CP="$CP:$jar"
done
echo "$CP" > bench_logs/.classpath.txt
echo "Classpath entries: $(echo $CP | tr ':' '\n' | wc -l)"

# 1. CCD baseline (no GPU needed, but keep same resources for fair comparison)
JID=$(sbatch $SLURM --job-name=bench_s7_ccd \
    --output=bench_logs/bench_s7_ccd_%j.out --error=bench_logs/bench_s7_ccd_%j.err \
    --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date)\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=ccd -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -cp \"\$(cat bench_logs/.classpath.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted CCD baseline: $JID"

# 2. Strategy7 with varying gpuBatchSize
for BS in "${BATCH_SIZES[@]}"; do
    JID=$(sbatch $SLURM --job-name=bench_s7_bs${BS} \
        --output=bench_logs/bench_s7_bs${BS}_%j.out --error=bench_logs/bench_s7_bs${BS}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date), Strategy7 gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy7 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')
    echo "Submitted Strategy7 gpuBatch=$BS: $JID"
done

echo ""
TOTAL=$((1 + ${#BATCH_SIZES[@]}))
echo "All $TOTAL jobs submitted (parallel). Monitor: squeue -u \$USER"
