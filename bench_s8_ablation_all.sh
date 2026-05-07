#!/bin/bash
# Submit all 5 S8 architecture benchmarks with dependency on round1
# 8 GPUs on fennario-06, 2 per task → 4 concurrent, SLURM will queue the 5th

cd /home/users/lz280/IdeaProjects/OSPREY3

SLURM="--partition=grisman --account=grisman --nodelist=fennario-06 --gres=gpu:a5000:2 --cpus-per-task=8 --mem=32G --time=02:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"
DEPEND="--dependency=afterany:11124526:11124527"

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
ENV_SETUP="export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH"

declare -A SUFFIX_MAP
SUFFIX_MAP[baseline]="model_subtree"
SUFFIX_MAP[transformer]="model_subtree_transformer"
SUFFIX_MAP[dual_enc]="model_subtree_dual_enc"
SUFFIX_MAP[multihead]="model_subtree_multihead"
SUFFIX_MAP[transfer]="model_subtree_transfer"

for ARCH in baseline transformer dual_enc multihead transfer; do
    SUFFIX="${SUFFIX_MAP[$ARCH]}"
    PSUBTREE="$BASE/protein/$SUFFIX/subtree_model.onnx"
    CSUBTREE="$BASE/complex/$SUFFIX/subtree_model.onnx"

    sbatch $SLURM $DEPEND --job-name=bench_s8_${ARCH} \
        --output=bench_logs/bench_s8_${ARCH}_%j.out --error=bench_logs/bench_s8_${ARCH}_%j.err \
        --wrap "$ENV_SETUP && echo \"Node: \$(hostname), Date: \$(date), S8 arch=$ARCH gpuBatch=$BS\" && $JAVA $JARGS -Dosprey.gnn.benchStrategy=strategy8 -Dosprey.gnn.confSpace=$CONFSPACE -Dosprey.gnn.numSeqs=$NUMSEQS -Dosprey.gnn.gpuBatchSize=$BS -Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.eval.proteinSubtreeModelPath=$PSUBTREE -Dosprey.gnn.eval.complexSubtreeModelPath=$CSUBTREE -Dosprey.gnn.trainingConfSpace=all20 -cp \"\$(cat bench_logs/.classpath_s8abl.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\""
done

echo "5 S8 jobs submitted (depend on 11124526,11124527)"
