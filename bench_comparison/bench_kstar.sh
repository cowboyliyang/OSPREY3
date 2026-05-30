#!/bin/bash
# ===================================================================
# K* baseline benchmark on dance_bench PDBs
# Usage: bash bench_kstar.sh [design_id or "all"]
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
DESIGN=${1:-d004}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx64g -Xms4g"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
SLURM="--partition=compsci --cpus-per-task=16 --mem=64G --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"

mkdir -p "$LOGDIR" "$OUTDIR/results"

# Build & classpath
./gradlew testClasses 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
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

    # Auto-detect chains: first chain = protein, rest = ligand (simplified)
    # For real runs, may need per-PDB chain config
    JID=$(sbatch $SLURM --job-name=kstar_${did} \
        --output=$LOGDIR/kstar_${did}_%j.out --error=$LOGDIR/kstar_${did}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=kstar \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results \
            -Dosprey.bench.numCPUs=16 \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    echo "Submitted K* $did ($pdb): $JID"
}

if [ "$DESIGN" = "all" ]; then
    for did in $(grep -v "^#" "$SPECS" | grep -v "^[[:space:]]*$" | cut -d',' -f1); do submit_one "$did"; done
else
    submit_one "$DESIGN"
fi
