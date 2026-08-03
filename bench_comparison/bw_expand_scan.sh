#!/bin/bash
# ===================================================================
# Branchwidth-expansion scan: for a given design, incrementally add
# interface-proximal residues (closest first, from rank_interface_candidates.py)
# to the flexible set, running PACK* dry-run branchwidth diagnostics at each
# step, until branchwidth clearly exceeds the omega=11 ceiling.
# Cutoff strategy: RESIDUAL_BUDGET at residualBudget=3.0 (per instruction).
# Records (k_added, branchwidth) so we can later pick the max k with
# branchwidth<=10 and the max k with branchwidth<=11.
#
# Usage: bash bw_expand_scan.sh <design_id>
# ===================================================================
set -uo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

DESIGN=$1
OUTDIR=/usr/xtmp/lz280/bench_comparison/bw_expand
SPECS=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
RESID_BUDGET=${RESID_BUDGET:-3.0}
MAX_ADD=${MAX_ADD:-80}
STOP_MARGIN=${STOP_MARGIN:-13}       # stop once bw exceeds this for STOP_STREAK consecutive steps
STOP_STREAK=${STOP_STREAK:-4}

mkdir -p "$OUTDIR"
RESULT="$OUTDIR/${DESIGN}_scan.tsv"
LOG="$OUTDIR/${DESIGN}_scan.log"
: > "$RESULT"
: > "$LOG"

line=$(grep "^${DESIGN}," "$SPECS")
if [ -z "$line" ]; then echo "Design $DESIGN not found in specs"; exit 1; fi
pdb=$(echo "$line" | cut -d',' -f2)
mutable=$(echo "$line" | cut -d, -f6)
flexible=$(echo "$line" | cut -d, -f7)
pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"

echo "[$DESIGN] mutable=$mutable flexible(orig)=$flexible" | tee -a "$LOG"

CAND="$OUTDIR/${DESIGN}_candidates.txt"
python3 bench_comparison/rank_interface_candidates.py "$pdbpath" "$mutable" "$flexible" "$CAND" 2>&1 | tee -a "$LOG"
NCAND=$(wc -l < "$CAND")
echo "[$DESIGN] $NCAND candidates ranked" | tee -a "$LOG"

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP=$(cat /usr/xtmp/lz280/bench_comparison/logs/.classpath_bench.txt 2>/dev/null)
if [ -z "$CP" ]; then
    echo "[$DESIGN] no cached classpath, building..." | tee -a "$LOG"
    ./gradlew testClasses 2>&1 | tail -5
    CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
    for jar in lib/*.jar; do CP="$CP:$jar"; done
    for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
fi

JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx48g -Xms4g"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

streak=0
k=0
echo -e "k\tbranchwidth\tadded_residues" > "$RESULT"

# k=0 baseline (original design, no additions) first
run_one() {
    local k=$1
    local extra_flex=$2
    local full_flex="$flexible"
    if [ -n "$extra_flex" ]; then full_flex="$flexible;$extra_flex"; fi

    local out
    out=$($JAVA $JARGS \
        -Dbranchdp.dp.dryRun=true \
        -Dbranchmarkstar.dp.cache=false \
        -Dosprey.bench.method=pac \
        -Dbranchmarkstar.usePAC=true \
        -Dbranchdp.rootSplit=legacy \
        -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET \
        -Dbranchdp.cutoff.residualBudget=$RESID_BUDGET \
        -Dosprey.bench.pdbPath="$pdbpath" \
        -Dosprey.bench.mutable="$mutable" \
        -Dosprey.bench.flexible="$full_flex" \
        -Dosprey.bench.method=pac \
        -Dosprey.bench.designId=$DESIGN \
        -Dosprey.bench.outputDir=/tmp/bw_expand_dummy_$DESIGN \
        -Dosprey.bench.numCPUs=8 \
        -cp "$CP" $MAIN 2>&1)

    echo "$out" >> "$LOG"
    # take max branchwidth reported across protein/ligand/complex blocks
    local bw
    bw=$(echo "$out" | grep -oP "Branchwidth=\K[0-9]+" | sort -n | tail -1)
    if [ -z "$bw" ]; then bw=-1; fi
    echo -e "${k}\t${bw}\t${extra_flex}" >> "$RESULT"
    echo "[$DESIGN] k=$k bw=$bw" | tee -a "$LOG"
    echo "$bw"
}

bw0=$(run_one 0 "")

added=""
for ((k=1; k<=MAX_ADD && k<=NCAND; k++)); do
    res=$(sed -n "${k}p" "$CAND" | cut -f1)
    if [ -z "$added" ]; then added="$res"; else added="$added;$res"; fi
    bw=$(run_one "$k" "$added")
    if [ "$bw" -gt "$STOP_MARGIN" ] 2>/dev/null; then
        streak=$((streak+1))
        if [ "$streak" -ge "$STOP_STREAK" ]; then
            echo "[$DESIGN] stopping early at k=$k (bw>$STOP_MARGIN for $STOP_STREAK consecutive steps)" | tee -a "$LOG"
            break
        fi
    else
        streak=0
    fi
done

echo "[$DESIGN] DONE" | tee -a "$LOG"
