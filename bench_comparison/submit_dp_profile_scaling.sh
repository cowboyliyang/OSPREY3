#!/bin/bash
# Submit BranchMARK* DP-only scaling profiles for the known-large 2xxm case.
#
# Defaults reproduce the completed PAC run's large final sequence:
#   design=2xxm, state=complex, seqIndex=19, rootSplit=legacy
#
# Env overrides:
#   DESIGN=2xxm STATE=complex SEQ_INDEX=19 REPS=3
#   THREADS_LIST="1 2 4 8 16 32"
#   PARTITION=compsci CPUS=128 MEM=256G TIME=08:00:00 NODES=1
#   ACCOUNT=grisman NODELIST=fennario-01,fennario-02 CONSTRAINT=a5000 EXCLUSIVE=false
#   JAVA_XMX=192g JAVA_XMS=4g ROOT_SPLIT=legacy

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=${OUTDIR:-/usr/xtmp/lz280/bench_comparison}
LOGDIR=${LOGDIR:-$OUTDIR/logs/dp_profile_2xxm}
SPECS=${SPECS:-$OUTDIR/design_specs_prepped.csv}
PDBDIR=${PDBDIR:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}

DESIGN=${DESIGN:-2xxm}
STATE=${STATE:-complex}
SEQ_INDEX=${SEQ_INDEX:-19}
MAX_MUT=${MAX_MUT:-1}
REPS=${REPS:-3}
THREADS_LIST=${THREADS_LIST:-"1 2 4 8 16 32"}

PARTITION=${PARTITION:-compsci}
CPUS=${CPUS:-128}
MEM=${MEM:-256G}
TIME=${TIME:-08:00:00}
NODES=${NODES:-1}
ACCOUNT=${ACCOUNT:-}
NODELIST=${NODELIST:-}
CONSTRAINT=${CONSTRAINT:-}
EXCLUSIVE=${EXCLUSIVE:-false}

ROOT_SPLIT=${ROOT_SPLIT:-legacy}
DP_CACHE=${DP_CACHE:-false}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-4g}
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

mkdir -p "$LOGDIR" "$OUTDIR/results"

line=$(grep "^${DESIGN}," "$SPECS" || true)
if [ -z "$line" ]; then
    echo "Design $DESIGN not found in $SPECS" >&2
    exit 1
fi

pdb=$(echo "$line" | cut -d',' -f2)
mutable=$(echo "$line" | cut -d',' -f6)
flexible=$(echo "$line" | cut -d',' -f7)
pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"

if [ ! -f "$pdbpath" ]; then
    echo "PDB not ready: $pdbpath" >&2
    exit 1
fi

CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do
    CP="$CP:$jar"
done
while IFS= read -r jar; do
    CP="$CP:$jar"
done < <(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u)
echo "$CP" > "$LOGDIR/.classpath_dp_profile.txt"

JARGS="--add-opens java.base/java.util=ALL-UNNAMED \
--add-opens java.base/java.lang=ALL-UNNAMED \
--add-opens java.base/java.lang.invoke=ALL-UNNAMED \
-Xmx$JAVA_XMX -Xms$JAVA_XMS -XX:-UseSuperWord"

echo "Submitting DP profile scaling:"
echo "  design=$DESIGN pdb=$pdb state=$STATE seqIndex=$SEQ_INDEX"
echo "  mutable=$mutable"
echo "  flexible=$flexible"
echo "  partition=$PARTITION nodes=$NODES cpus=$CPUS mem=$MEM time=$TIME rootSplit=$ROOT_SPLIT"
echo "  account=${ACCOUNT:-none} nodelist=${NODELIST:-none} constraint=${CONSTRAINT:-none} exclusive=$EXCLUSIVE"
echo "  threads=$THREADS_LIST reps=$REPS"
echo "  logs=$LOGDIR"

SBATCH_FLAGS=(--parsable
    --partition="$PARTITION"
    --time="$TIME"
    --nodes="$NODES"
    --cpus-per-task="$CPUS"
    --job-name=placeholder)
if [ "$EXCLUSIVE" = "true" ]; then
    SBATCH_FLAGS+=(--exclusive --mem=0)
else
    SBATCH_FLAGS+=(--mem="$MEM")
fi
if [ -n "$ACCOUNT" ]; then
    SBATCH_FLAGS+=(--account="$ACCOUNT")
fi
if [ -n "$NODELIST" ]; then
    SBATCH_FLAGS+=(--nodelist="$NODELIST")
fi
if [ -n "$CONSTRAINT" ]; then
    SBATCH_FLAGS+=(--constraint="$CONSTRAINT")
fi

for threads in $THREADS_LIST; do
    for rep in $(seq 1 "$REPS"); do
        flags=("${SBATCH_FLAGS[@]}")
        for i in "${!flags[@]}"; do
            if [ "${flags[$i]}" = "--job-name=placeholder" ]; then
                flags[$i]="--job-name=dp2xxm_t${threads}_r${rep}"
            fi
        done
        jid=$(sbatch "${flags[@]}" \
            --output="$LOGDIR/dp_${DESIGN}_${STATE}_seq${SEQ_INDEX}_t${threads}_r${rep}_%j.out" \
            --error="$LOGDIR/dp_${DESIGN}_${STATE}_seq${SEQ_INDEX}_t${threads}_r${rep}_%j.err" \
            --wrap "$JAVA $JARGS \
                -Dosprey.bench.pdbPath=$pdbpath \
                -Dosprey.bench.mutable='$mutable' \
                -Dosprey.bench.flexible='$flexible' \
                -Dosprey.bench.method=dp_profile \
                -Dosprey.bench.designId=$DESIGN \
                -Dosprey.bench.outputDir=$OUTDIR/results \
                -Dosprey.bench.numCPUs=$CPUS \
                -Dosprey.dpProfile.state=$STATE \
                -Dosprey.dpProfile.seqIndex=$SEQ_INDEX \
                -Dosprey.dpProfile.maxMut=$MAX_MUT \
                -Dbranchmarkstar.rootSplit=$ROOT_SPLIT \
                -Dbranchmarkstar.dp.cache=$DP_CACHE \
                -Dbranchmarkstar.dp.parallel=true \
                -Dbranchmarkstar.dp.parallel.minMStates=1 \
                -Dbranchmarkstar.dp.parallel.threads=$threads \
                -Dbranchmarkstar.dp.progress=true \
                -cp \"\$(cat $LOGDIR/.classpath_dp_profile.txt)\" $MAIN 2>&1")
        echo "submitted threads=$threads rep=$rep job=$jid"
    done
done
