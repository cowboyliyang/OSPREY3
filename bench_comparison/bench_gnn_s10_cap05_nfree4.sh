#!/bin/bash
# ===================================================================
# MARK* + GNN S9 benchmark on dance_bench PDBs
# Requires trained GNN models. Run train_gnn.sh first.
# Usage: bash bench_gnn_s9.sh [design_id or "all"]
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
GNNDIR=$OUTDIR/gnn_models
DESIGN=${1:-d004}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ORT_CUDA_LIB=/home/users/lz280/IdeaProjects/OSPREY3/lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx64g -Xms4g -Djava.library.path=$ORT_CUDA_LIB:$CUDA_LIB"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
SLURM="--partition=grisman --account=grisman --gres=gpu:a5000:1 --cpus-per-task=8 --mem=64G --time=7-00:00:00 --nodelist=fennario-01,fennario-02,fennario-03 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

mkdir -p "$LOGDIR" "$OUTDIR/results/gnn_s10_cap05_nfree4"

./gradlew testClasses 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_bench.txt"

submit_one() {
    local did=$1
    local line=$(grep "^${did}," "$SPECS")
    if [ -z "$line" ]; then echo "Design $did not found"; return; fi

    local pdb=$(echo "$line" | cut -d',' -f2)
    local mutable=$(echo "$line" | cut -d, -f6)
    local flexible=$(echo "$line" | cut -d, -f7)
    local pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"

    if [ ! -f "$pdbpath" ]; then echo "PDB not ready: $pdbpath"; return; fi

    # Pass only models that actually exist; Java side handles per-confspace fallback.
    local pModel="$GNNDIR/$did/protein/model/gnn_model.onnx"
    local cModel="$GNNDIR/$did/complex/model/gnn_model.onnx"
    local lModel="$GNNDIR/$did/ligand/model/gnn_model.onnx"
    local pSub="$GNNDIR/$did/protein/model_subtree/subtree_model.onnx"
    local cSub="$GNNDIR/$did/complex/model_subtree/subtree_model.onnx"
    local lSub="$GNNDIR/$did/ligand/model_subtree/subtree_model.onnx"

    local GNN_PROPS=""
    [ -f "$pModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.proteinModelPath=$pModel"
    [ -f "$cModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.complexModelPath=$cModel"
    [ -f "$lModel" ] && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.ligandModelPath=$lModel"
    [ -f "$pSub" ]   && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.proteinSubtreeModelPath=$pSub"
    [ -f "$cSub" ]   && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.complexSubtreeModelPath=$cSub"
    [ -f "$lSub" ]   && GNN_PROPS="$GNN_PROPS -Dosprey.gnn.eval.ligandSubtreeModelPath=$lSub"

    if [ -z "$GNN_PROPS" ]; then
        echo "No GNN models for $did — skipping (would be identical to markstar run)"
        return
    fi

    JID=$(sbatch $SLURM --job-name=s10n4_${did} \
        --output=$LOGDIR/gnn_s10_cap05_nfree4_${did}_%j.out --error=$LOGDIR/gnn_s10_cap05_nfree4_${did}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && export LD_LIBRARY_PATH=$ORT_CUDA_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH && $JAVA $JARGS \
            -Dmarkstar.fullParallelFromStart=true \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=gnn_s10 \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results/gnn_s10_cap05_nfree4 \
            -Dosprey.bench.numCPUs=8 \
            $GNN_PROPS \
            -Dosprey.gnn.gpuBatchSize=1000 \
            -Dosprey.gnn.s10.tauLeaf=0.3 \
            -Dosprey.gnn.s10.tauSubtree=0.5 \
            -Dosprey.gnn.s10.leafCpQ=0.3 \
            -Dosprey.gnn.s10.subtreeCpQ=0.5 \
            -Dosprey.gnn.s10.subtreeFinalize=true \
            -Dosprey.gnn.s10.maxBernstein=10000 \
            -Dosprey.gnn.s10.maxSubtreeFinalizeFree=4 \
            -Dosprey.gnn.s10.maxSubtreeFinalizeCpQ=0.5 \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    echo "Submitted GNN S10 $did ($pdb): $JID"
}

if [ "$DESIGN" = "all" ]; then
    for did in $(grep -v "^#" "$SPECS" | grep -v "^[[:space:]]*$" | cut -d',' -f1); do submit_one "$did"; done
else
    submit_one "$DESIGN"
fi
