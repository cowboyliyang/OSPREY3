#!/bin/bash
# ===================================================================
# BranchMARK* PAC-mode benchmark on dance_bench PDBs
#
# PAC = Probably Approximately Correct partition-function estimation via
# Rao-Blackwellized importance sampling on the BranchMARK* tree.
# It is CPU-bound in both DP-proposal sampling (Phase 1) and CCD minimization
# (Phase 2), so we run on the highest-core nodes available. CPU node survey:
#   grisman:  fennario-01..06 = 104 CPUs (+8x A5000)   <-- used here
#             jerry1..7 / grisman-37/40 = only 48 CPUs
#   compsci:  compsci-cluster-fitz-35..44 = 128 CPUs   (alt, set PARTITION=compsci)
#
# Runs on BOTH partitions, each job maxing out its node's cores ("拉满"):
#   grisman:  fennario-01..06 = 104 CPUs (+8x A5000, --constraint=a5000)
#   compsci:  compsci-cluster-fitz-35..44 = 128 CPUs (pinned via --cpus-per-task=128)
# For "all", designs are round-robined across the two partitions for max throughput.
#
# Usage: bash bench_pac.sh [design_id or "all"]
#   env overrides: PAC_SAMPLES=500  PAC_CONFIDENCE=0.05
#                  PAC_ADAPTIVE=true PAC_MIN_SAMPLES=100 PAC_BATCH_SIZE=200
#                  PAC_SAMPLING_BATCHED=true PAC_SAMPLING_THREADS= PAC_SAMPLING_LARGE_LAMBDA=65536
#                  ROOT_SPLIT=memory
#                  GRISMAN_CPUS=104  COMPSCI_CPUS=128
#                  GRISMAN_MEM=128G   COMPSCI_MEM=256G
#                  JAVA_XMX=192g      JAVA_XMS=4g
#                  DP_CACHE=false     EXTRA_JVM_ARGS="..."
#                  GRISMAN_FULL_NODE=false GRISMAN_EXCLUDE="" GRISMAN_CONSTRAINT=a5000
#                  (set GRISMAN_CONSTRAINT= to allow any grisman node)
#                  TARGET=both|grisman|compsci   (single-design default: grisman)
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
DESIGN=${1:-2q1e}

# --- tunables ---
PAC_SAMPLES=${PAC_SAMPLES:-500}        # branchmarkstar.pac.samples
PAC_CONFIDENCE=${PAC_CONFIDENCE:-0.05} # delta; 0.05 => 95% confidence
PAC_ADAPTIVE=${PAC_ADAPTIVE:-true}
PAC_MIN_SAMPLES=${PAC_MIN_SAMPLES:-100}
PAC_MAX_SAMPLES=${PAC_MAX_SAMPLES:-$PAC_SAMPLES}
PAC_BATCH_SIZE=${PAC_BATCH_SIZE:-200}
PAC_TARGET_EPSILON=${PAC_TARGET_EPSILON:-}
PAC_SAMPLING_BATCHED=${PAC_SAMPLING_BATCHED:-true}
PAC_SAMPLING_THREADS=${PAC_SAMPLING_THREADS:-}
PAC_SAMPLING_LARGE_LAMBDA=${PAC_SAMPLING_LARGE_LAMBDA:-65536}
GRISMAN_CPUS=${GRISMAN_CPUS:-104}      # max out fennario
COMPSCI_CPUS=${COMPSCI_CPUS:-128}      # max out fitz
GRISMAN_MEM=${GRISMAN_MEM:-128G}
COMPSCI_MEM=${COMPSCI_MEM:-256G}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-4g}
DP_CACHE=${DP_CACHE:-false}
ROOT_SPLIT=${ROOT_SPLIT:-memory}
GRISMAN_FULL_NODE=${GRISMAN_FULL_NODE:-false}
GRISMAN_EXCLUDE=${GRISMAN_EXCLUDE:-}
GRISMAN_CONSTRAINT=${GRISMAN_CONSTRAINT-a5000}
TARGET=${TARGET:-both}                 # both | grisman | compsci

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord"
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx$JAVA_XMX -Xms$JAVA_XMS $EXTRA_JVM_ARGS"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

# Per-partition SLURM flags. grisman pins the 104-CPU a5000 (fennario) nodes;
# compsci asking for 128 cpus-per-task already restricts to the fitz nodes.
slurm_flags() {
    local part=$1 cpus=$2 mem=$3
    local f="--partition=$part --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"
    if [ "$part" = "grisman" ]; then
        f="$f --account=grisman"
        if [ "$GRISMAN_FULL_NODE" = "true" ]; then
            f="$f --exclusive --mem=0"
        else
            f="$f --cpus-per-task=$cpus --mem=$mem"
        fi
        if [ -n "$GRISMAN_CONSTRAINT" ]; then f="$f --constraint=$GRISMAN_CONSTRAINT"; fi
        if [ -n "$GRISMAN_EXCLUDE" ]; then f="$f --exclude=$GRISMAN_EXCLUDE"; fi
    else
        f="$f --cpus-per-task=$cpus --mem=$mem"
    fi
    echo "$f"
}

mkdir -p "$LOGDIR" "$OUTDIR/results"

./gradlew testClasses 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_bench.txt"

submit_one() {
    local did=$1 part=$2 cpus=$3
    local line=$(grep "^${did}," "$SPECS")
    if [ -z "$line" ]; then echo "Design $did not found"; return; fi

    local pdb=$(echo "$line" | cut -d',' -f2)
    local mutable=$(echo "$line" | cut -d, -f6)
    local flexible=$(echo "$line" | cut -d, -f7)
    local pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"

    if [ ! -f "$pdbpath" ]; then echo "PDB not ready: $pdbpath"; return; fi

    local mem=$GRISMAN_MEM
    if [ "$part" = "compsci" ]; then mem=$COMPSCI_MEM; fi

    JID=$(sbatch $(slurm_flags "$part" "$cpus" "$mem") --job-name=pac_${did} \
        --output=$LOGDIR/pac_${did}_%j.out --error=$LOGDIR/pac_${did}_%j.err \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-$cpus}; $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=pac \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results \
            -Dosprey.bench.numCPUs=\$RUN_CPUS \
            -Dbranchmarkstar.usePAC=true \
            -Dbranchmarkstar.rootSplit=$ROOT_SPLIT \
            -Dbranchmarkstar.dp.cache=$DP_CACHE \
            -Dbranchmarkstar.pac.samples=$PAC_SAMPLES \
            -Dbranchmarkstar.pac.confidence=$PAC_CONFIDENCE \
            -Dbranchmarkstar.pac.adaptive=$PAC_ADAPTIVE \
            -Dbranchmarkstar.pac.minSamples=$PAC_MIN_SAMPLES \
            -Dbranchmarkstar.pac.maxSamples=$PAC_MAX_SAMPLES \
            -Dbranchmarkstar.pac.batchSize=$PAC_BATCH_SIZE \
            -Dbranchmarkstar.pac.targetEpsilon=$PAC_TARGET_EPSILON \
            -Dbranchmarkstar.pac.sampling.batched=$PAC_SAMPLING_BATCHED \
            -Dbranchmarkstar.pac.sampling.threads=$PAC_SAMPLING_THREADS \
            -Dbranchmarkstar.pac.sampling.largeLambdaThreshold=$PAC_SAMPLING_LARGE_LAMBDA \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    if [ "$part" = "grisman" ] && [ "$GRISMAN_FULL_NODE" = "true" ]; then
        echo "Submitted PAC $did ($pdb) [$part, full-node, exclude=${GRISMAN_EXCLUDE:-none}, constraint=${GRISMAN_CONSTRAINT:-none}, Xmx=${JAVA_XMX}]: $JID"
    else
        echo "Submitted PAC $did ($pdb) [$part, ${cpus} CPUs, ${mem}, Xmx=${JAVA_XMX}]: $JID"
    fi
}

if [ "$DESIGN" = "all" ]; then
    i=0
    for did in $(grep -v "^#" "$SPECS" | grep -v "^[[:space:]]*$" | cut -d',' -f1); do
        case "$TARGET" in
            grisman) submit_one "$did" grisman "$GRISMAN_CPUS" ;;
            compsci) submit_one "$did" compsci "$COMPSCI_CPUS" ;;
            both)    if [ $((i % 2)) -eq 0 ]; then submit_one "$did" grisman "$GRISMAN_CPUS";
                     else submit_one "$did" compsci "$COMPSCI_CPUS"; fi ;;
        esac
        i=$((i + 1))
    done
else
    # single design: default to grisman unless TARGET=compsci
    if [ "$TARGET" = "compsci" ]; then submit_one "$DESIGN" compsci "$COMPSCI_CPUS";
    else submit_one "$DESIGN" grisman "$GRISMAN_CPUS"; fi
fi
