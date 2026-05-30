#!/bin/bash
# ===================================================================
# GNN training for a specific design
# Step 1: Export training data via OSPREY (Java)
# Step 2: Train leaf + subtree GNN (Python)
#
# Usage: bash train_gnn.sh <pdb_id>    e.g. bash train_gnn.sh 4hem
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
GNNDIR=$OUTDIR/gnn_models
DESIGN=${1:?Usage: train_gnn.sh <pdb_id>}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx64g -Xms4g"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
SLURM="--partition=grisman --account=grisman --gres=gpu:a5000:1 --cpus-per-task=8 --mem=64G --time=12:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

mkdir -p "$LOGDIR"

# Read design spec (CSV format: pdb_id,pdb,desc,num_seq,num_flex_mut,mutable,flexible,...)
line=$(grep "^${DESIGN}," "$SPECS")
if [ -z "$line" ]; then echo "Design $DESIGN not found in $SPECS"; exit 1; fi

pdb=$(echo "$line" | cut -d',' -f2)
mutable=$(echo "$line" | cut -d',' -f6)
flexible=$(echo "$line" | cut -d',' -f7)
pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.pdb"
if [ ! -f "$pdbpath" ]; then echo "PDB not ready: $pdbpath"; exit 1; fi

./gradlew testClasses 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_bench.txt"

DATA_DIR="$GNNDIR/$DESIGN"
mkdir -p "$DATA_DIR"

# Submit combined export + train job
JID=$(sbatch $SLURM --job-name=train_${DESIGN} \
    --output=$LOGDIR/train_${DESIGN}_%j.out --error=$LOGDIR/train_${DESIGN}_%j.err \
    --wrap "
cd /home/users/lz280/IdeaProjects/OSPREY3

echo '=== Step 1: Export GNN training data ==='
echo 'Design: $DESIGN, PDB: $pdb'
echo 'Mutable: $mutable, Flexible: $flexible'

# Export training data (Java — generic, any PDB)
$JAVA $JARGS \
    -Dosprey.bench.pdbPath=$pdbpath \
    -Dosprey.bench.mutable='$mutable' \
    -Dosprey.bench.flexible='$flexible' \
    -Dosprey.bench.method=export_gnn \
    -Dosprey.bench.designId=$DESIGN \
    -Dosprey.bench.outputDir=$OUTDIR/results \
    -Dosprey.gnn.outputDir=$DATA_DIR \
    -Dosprey.gnn.numSamples=200000 \
    -Dosprey.bench.numCPUs=8 \
    -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1

echo ''
echo '=== Step 2: Train GNN models ==='
eval \"\$(/home/users/lz280/miniconda3/bin/conda shell.bash hook)\"
conda activate confdiff

for space in protein complex; do
    data=$DATA_DIR/\$space
    if [ ! -f \"\$data/confs.csv\" ]; then
        echo \"SKIP \$space: no confs.csv\"
        continue
    fi
    echo \"--- Training leaf GNN: \$space ---\"
    python gnn/train.py --data_dir \"\$data\" --epochs 200 --batch_size 512

    echo \"--- Training subtree GNN: \$space ---\"
    python gnn/train_subtree.py --data_dir \"\$data\" --epochs 200 --batch_size 256
done

echo ''
echo '=== Training complete ==='
ls -la $DATA_DIR/*/model*/*.onnx 2>/dev/null || echo 'WARNING: no ONNX models found'
" | awk '{print $4}')
echo "Submitted GNN training for $DESIGN ($pdb): $JID"
