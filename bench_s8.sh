#!/bin/bash
# Submit GNN benchmark: CCD baseline vs Strategy8 (S7 + subtree GNN for internal nodes)
# Uses A5000 on grisman partition.
# Usage: bash bench_s8.sh [gpuBatchSize1 gpuBatchSize2 ...]
# Default batch sizes: 500, 1000, 2000

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

SLURM="--partition=grisman --account=grisman --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=02:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

# Per-conf GNN models (leaf nodes)
PMODEL="gnn_data/2RL0_all20_4pos_merged/protein/model/gnn_model.onnx"
CMODEL="gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_model.onnx"

# Subtree GNN models (internal nodes)
PSUBTREE="gnn_data/2RL0_all20_4pos_merged/protein/model_subtree/subtree_model.onnx"
CSUBTREE="gnn_data/2RL0_all20_4pos_merged/complex/model_subtree/subtree_model.onnx"

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
echo "$CP" > bench_logs/.classpath_s8.txt
echo "Classpath entries: $(echo $CP | tr ':' '\n' | wc -l)"

# Check that subtree models exist
for f in "$PSUBTREE" "$CSUBTREE"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Subtree model not found: $f"
        echo "  Run subtree GNN training first, then re-export ONNX."
        exit 1
    fi
done

# 1. CCD baseline
JID=$(sbatch $SLURM --job-name=bench_s8_ccd \
    --output=bench_logs/bench_s8_ccd_%j.out --error=bench_logs/bench_s8_ccd_%j.err \
    --wrap "export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date)\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=ccd -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -cp \"\$(cat bench_logs/.classpath_s8.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted CCD baseline: $JID"

# 2. Strategy7 baseline (for comparison)
for BS in "${BATCH_SIZES[@]}"; do
    JID=$(sbatch $SLURM --job-name=bench_s7_bs${BS} \
        --output=bench_logs/bench_s8_s7ref_bs${BS}_%j.out --error=bench_logs/bench_s8_s7ref_bs${BS}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date), Strategy7 gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy7 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')
    echo "Submitted Strategy7 gpuBatch=$BS: $JID"
done

# 3. Strategy8 with varying gpuBatchSize
for BS in "${BATCH_SIZES[@]}"; do
    JID=$(sbatch $SLURM --job-name=bench_s8_bs${BS} \
        --output=bench_logs/bench_s8_bs${BS}_%j.out --error=bench_logs/bench_s8_bs${BS}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date), Strategy8 gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy8 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.eval.proteinSubtreeModelPath=$PSUBTREE -Dosprey.gnn.eval.complexSubtreeModelPath=$CSUBTREE -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')
    echo "Submitted Strategy8 gpuBatch=$BS: $JID"
done

# 4. Strategy9 with varying gpuBatchSize (logZ residual subtree bound oracle)
for BS in "${BATCH_SIZES[@]}"; do
    JID=$(sbatch $SLURM --job-name=bench_s9_bs${BS} \
        --output=bench_logs/bench_s9_bs${BS}_%j.out --error=bench_logs/bench_s9_bs${BS}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date), Strategy9 gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy9 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.eval.proteinSubtreeModelPath=$PSUBTREE -Dosprey.gnn.eval.complexSubtreeModelPath=$CSUBTREE -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')
    echo "Submitted Strategy9 gpuBatch=$BS: $JID"
done

TOTAL=$((1 + 3 * ${#BATCH_SIZES[@]}))
echo "All $TOTAL jobs submitted. Monitor: squeue -u \$USER"
