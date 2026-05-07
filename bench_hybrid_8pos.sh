#!/bin/bash
# Submit 4 hybrid 8-pos benchmarks (4 mutable highrot + 4 flexible):
#   1. markstar_ccd    — MARK* + CCD
#   2. markstar_s8     — MARK* + GNN Strategy 8
#   3. branch_ccd      — BranchMARK* + CCD
#   4. branch_s8       — BranchMARK* + GNN Strategy 8
#
# Uses A5000 GPUs on grisman partition (fennario nodes).
# GNN models: trained on full 8-pos highrot, with RC mapping for hybrid confspace.
#
# Usage: bash bench_hybrid_8pos.sh

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

# SLURM settings — A5000 on grisman
# GNN methods need 4 GPUs (protein leaf/subtree + complex leaf/subtree)
SLURM_GPU="--partition=grisman --account=grisman --gres=gpu:titan_v:4 --cpus-per-task=16 --mem=64G --time=12:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL --exclude=jerry1"
SLURM_CCD="--partition=grisman --account=grisman --cpus-per-task=16 --mem=64G --time=24:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL --exclude=jerry1"

# GNN models (8pos highrot)
PMODEL="gnn_data/2RL0_highrot_8pos/protein/model/gnn_model.onnx"
CMODEL="gnn_data/2RL0_highrot_8pos/complex/model/gnn_model.onnx"
PSUBTREE="gnn_data/2RL0_highrot_8pos/protein/model_subtree/subtree_model.onnx"
CSUBTREE="gnn_data/2RL0_highrot_8pos/complex/model_subtree/subtree_model.onnx"

EPSILON=0.683
CPUS=16
GPUBATCH=1000
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx60g -Xms4g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.TestBranchMARKStar

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
echo "$CP" > bench_logs/.classpath_hybrid8.txt
echo "Classpath entries: $(echo $CP | tr ':' '\n' | wc -l)"

# Check models exist
for f in "$PMODEL" "$CMODEL" "$PSUBTREE" "$CSUBTREE"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Model not found: $f"
        exit 1
    fi
done

# Common java properties
HYBRID_PROPS="-Dosprey.hybrid.epsilon=$EPSILON -Dosprey.hybrid.numCPUs=$CPUS"
GNN_PROPS="-Dosprey.gnn.eval.proteinModelPath=$PMODEL -Dosprey.gnn.eval.complexModelPath=$CMODEL -Dosprey.gnn.eval.proteinSubtreeModelPath=$PSUBTREE -Dosprey.gnn.eval.complexSubtreeModelPath=$CSUBTREE -Dosprey.gnn.gpuBatchSize=$GPUBATCH"

# 1. MARK* + CCD
JID1=$(sbatch $SLURM_CCD --job-name=hybrid_markstar_ccd \
    --output=bench_logs/hybrid_markstar_ccd_%j.out --error=bench_logs/hybrid_markstar_ccd_%j.err \
    --wrap "export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && echo \"Node: \$(hostname), Date: \$(date), Method: markstar_ccd\" && $JAVA $JARGS $HYBRID_PROPS -Dosprey.hybrid.method=markstar_ccd -cp \"\$(cat bench_logs/.classpath_hybrid8.txt)\" $MAIN hybrid_bench 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted markstar_ccd: $JID1"

# 2. MARK* + GNN S8
JID2=$(sbatch $SLURM_GPU --job-name=hybrid_markstar_s8 \
    --output=bench_logs/hybrid_markstar_s8_%j.out --error=bench_logs/hybrid_markstar_s8_%j.err \
    --wrap "export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && cd /home/users/lz280/IdeaProjects/OSPREY3 && echo \"Node: \$(hostname), Date: \$(date), Method: markstar_s8\" && $JAVA $JARGS $HYBRID_PROPS $GNN_PROPS -Dosprey.hybrid.method=markstar_s8 -cp \"\$(cat bench_logs/.classpath_hybrid8.txt)\" $MAIN hybrid_bench 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted markstar_s8: $JID2"

# 3. BranchMARK* + CCD
JID3=$(sbatch $SLURM_CCD --job-name=hybrid_branch_ccd \
    --output=bench_logs/hybrid_branch_ccd_%j.out --error=bench_logs/hybrid_branch_ccd_%j.err \
    --wrap "export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && cd /home/users/lz280/IdeaProjects/OSPREY3 && echo \"Node: \$(hostname), Date: \$(date), Method: branch_ccd\" && $JAVA $JARGS $HYBRID_PROPS -Dosprey.hybrid.method=branch_ccd -cp \"\$(cat bench_logs/.classpath_hybrid8.txt)\" $MAIN hybrid_bench 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted branch_ccd: $JID3"

# 4. BranchMARK* + GNN S8
JID4=$(sbatch $SLURM_GPU --job-name=hybrid_branch_s8 \
    --output=bench_logs/hybrid_branch_s8_%j.out --error=bench_logs/hybrid_branch_s8_%j.err \
    --wrap "export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && cd /home/users/lz280/IdeaProjects/OSPREY3 && echo \"Node: \$(hostname), Date: \$(date), Method: branch_s8\" && $JAVA $JARGS $HYBRID_PROPS $GNN_PROPS -Dosprey.hybrid.method=branch_s8 -cp \"\$(cat bench_logs/.classpath_hybrid8.txt)\" $MAIN hybrid_bench 2>&1 && echo \"Done: \$(date)\"" \
    | awk '{print $4}')
echo "Submitted branch_s8: $JID4"

echo ""
echo "All 4 jobs submitted. Monitor: squeue -u \$USER"
echo "  markstar_ccd: $JID1"
echo "  markstar_s8:  $JID2"
echo "  branch_ccd:   $JID3"
echo "  branch_s8:    $JID4"
