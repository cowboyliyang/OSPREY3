#!/bin/bash
# ===================================================================
# MARK* + GNN S11 leaf-only benchmark.
#
# S11 settings here:
#   - leaf GNN replaces CCD online
#   - subtree navigator is hard-disabled
#   - audit leaf CSVs are emitted in full mode
#
# Usage:
#   bash bench_comparison/bench_gnn_s11_leafonly.sh 3gxu
#   bash bench_comparison/bench_gnn_s11_leafonly.sh 8
#   bash bench_comparison/bench_gnn_s11_leafonly.sh all
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
GNNDIR=$OUTDIR/gnn_models
RESULT_DIR=$OUTDIR/results/gnn_s11_leafonly
AUDIT_DIR=$OUTDIR/audit_leaves/gnn_s11_leafonly
DESIGN=${1:-3gxu}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx64g -Xms4g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
SLURM="--partition=grisman --account=grisman --gres=gpu:a5000:1 --cpus-per-task=8 --mem=64G --time=7-00:00:00 --nodelist=fennario-01,fennario-02,fennario-03 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

mkdir -p "$LOGDIR" "$RESULT_DIR" "$AUDIT_DIR"

./gradlew testClasses 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_bench.txt"

submit_one() {
    local did=$1
    local line
    line=$(grep "^${did}," "$SPECS" || true)
    if [ -z "$line" ]; then echo "Design $did not found"; return; fi

    local pdb mutable flexible pdbpath
    pdb=$(echo "$line" | cut -d',' -f2)
    mutable=$(echo "$line" | cut -d, -f6)
    flexible=$(echo "$line" | cut -d, -f7)
    pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"

    if [ ! -f "$pdbpath" ]; then echo "PDB not ready: $pdbpath"; return; fi

    local pModel="$GNNDIR/$did/protein/model/gnn_model.onnx"
    local cModel="$GNNDIR/$did/complex/model/gnn_model.onnx"
    local lModel="$GNNDIR/$did/ligand/model/gnn_model.onnx"

    local GNN_PROPS=""
    [ -f "$pModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.proteinModelPath=$pModel"
    [ -f "$cModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.complexModelPath=$cModel"
    [ -f "$lModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.ligandModelPath=$lModel"

    if [ -z "$GNN_PROPS" ]; then
        echo "No leaf GNN models for $did -- skipping"
        return
    fi

    local JID
    JID=$(sbatch $SLURM --job-name=s11leaf_${did} \
        --output=$LOGDIR/gnn_s11_leafonly_${did}_%j.out --error=$LOGDIR/gnn_s11_leafonly_${did}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && $JAVA $JARGS \
            -Dmarkstar.fullParallelFromStart=true \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=gnn_s11 \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$RESULT_DIR \
            -Dosprey.bench.numCPUs=8 \
            $GNN_PROPS \
            -Dosprey.gnn.gpuBatchSize=1000 \
            -Dosprey.gnn.s11.subtreeNavigator=false \
            -Dosprey.gnn.s11.landscapeMode=off \
            -Dosprey.gnn.s11.auditMode=full \
            -Dosprey.gnn.s11.auditDir=$AUDIT_DIR \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    echo "Submitted GNN S11 leaf-only $did ($pdb): $JID"
}

if [ "$DESIGN" = "all" ]; then
    for did in $(grep -v "^#" "$SPECS" | grep -v "^[[:space:]]*$" | cut -d',' -f1); do submit_one "$did"; done
elif [ "$DESIGN" = "8" ]; then
    for did in 2p4a 2q1e 2rf9 3bu8 3gxu 3u7y 4pxf 4wem; do submit_one "$did"; done
else
    submit_one "$DESIGN"
fi
