#!/bin/bash
# Benchmark Strategy8 with different subtree GNN architectures on fennario-06
# 7 tasks total, 2x A5000 each, fennario-06 has 8 GPUs → 2 rounds
# Round 1 (4 jobs): CCD, S7-ref, S8-baseline, S8-transformer
# Round 2 (3 jobs): S8-dual_enc, S8-multihead, S8-transfer

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

SLURM="--partition=grisman --account=grisman --nodelist=fennario-06 --gres=gpu:a5000:2 --cpus-per-task=8 --mem=32G --time=02:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

# Leaf-node GNN models
PMODEL="gnn_data/2RL0_all20_4pos_merged/protein/model/gnn_model.onnx"
CMODEL="gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_model.onnx"

BASE="gnn_data/2RL0_all20_4pos_merged"
NUMSEQS=5
CONFSPACE=highrot
BS=1000
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx32g -Xms2g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.RunBenchmark

# Build
echo "Building..."
./gradlew testClasses 2>&1 | tail -3

# Classpath
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do
    CP="$CP:$jar"
done
echo "$CP" > bench_logs/.classpath_s8abl.txt
echo "Classpath entries: $(echo $CP | tr ':' '\n' | wc -l)"

ENV_SETUP="export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH"

submit_s8() {
    local ARCH=$1
    local SUFFIX=$2  # model_subtree or model_subtree_transformer etc.
    local PSUBTREE="$BASE/protein/$SUFFIX/subtree_model.onnx"
    local CSUBTREE="$BASE/complex/$SUFFIX/subtree_model.onnx"

    JID=$(sbatch $SLURM --job-name=bench_s8_${ARCH} \
        --output=bench_logs/bench_s8_${ARCH}_%j.out --error=bench_logs/bench_s8_${ARCH}_%j.err \
        --wrap "$ENV_SETUP && echo \"Node: \$(hostname), Date: \$(date), S8 arch=$ARCH gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy8 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.eval.proteinSubtreeModelPath=$PSUBTREE -Dosprey.gnn.eval.complexSubtreeModelPath=$CSUBTREE -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8abl.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')
    echo "  S8-$ARCH: $JID"
    echo "$JID"
}

# ============ Round 1: 4 jobs ============
echo ""
echo "=== Round 1 (4 jobs) ==="

# CCD baseline
CCD_JID=$(sbatch $SLURM --job-name=bench_s8abl_ccd \
    --output=bench_logs/bench_s8abl_ccd_%j.out --error=bench_logs/bench_s8abl_ccd_%j.err \
    --wrap "$ENV_SETUP && echo \"Node: \$(hostname), Date: \$(date)\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=ccd -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -cp \"\$(cat bench_logs/.classpath_s8abl.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "  CCD baseline: $CCD_JID"

# S7 reference
S7_JID=$(sbatch $SLURM --job-name=bench_s8abl_s7 \
    --output=bench_logs/bench_s8abl_s7_%j.out --error=bench_logs/bench_s8abl_s7_%j.err \
    --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && $ENV_SETUP && echo \"Node: \$(hostname), Date: \$(date), Strategy7 gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy7 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8abl.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "  S7 reference: $S7_JID"

# S8-baseline and S8-transformer
R1_J3=$(submit_s8 baseline model_subtree)
R1_J4=$(submit_s8 transformer model_subtree_transformer)

# Collect round 1 job IDs
R1_JIDS="$CCD_JID:$S7_JID:$R1_J3:$R1_J4"

# ============ Round 2: 3 jobs (after round 1 finishes) ============
echo ""
echo "=== Round 2 (3 jobs, after round 1) ==="

DEPEND="--dependency=afterany:${CCD_JID}:${S7_JID}:${R1_J3}:${R1_J4}"

# Override SLURM to add dependency
SLURM2="$SLURM $DEPEND"

# Temporarily swap SLURM for submit_s8
SLURM_ORIG="$SLURM"
SLURM="$SLURM2"
submit_s8 dual_enc model_subtree_dual_enc
submit_s8 multihead model_subtree_multihead
submit_s8 transfer model_subtree_transfer
SLURM="$SLURM_ORIG"

echo ""
echo "All 7 jobs submitted. Round 2 starts after round 1 finishes."
echo "Monitor: squeue -u \$USER"
