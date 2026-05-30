#!/bin/bash
# ===================================================================
# BranchMARK* PAC-mode benchmark on dance_bench PDBs
#
# PAC = Probably Approximately Correct partition-function estimation via
# Rao-Blackwellized importance sampling on the BranchMARK* tree.
# It is CPU-bound (parallel CCD sampling), so we run on the highest-core
# nodes available. CPU node survey:
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
#                  GRISMAN_CPUS=104  COMPSCI_CPUS=128
#                  GRISMAN_MEM=128G   COMPSCI_MEM=256G
#                  JAVA_XMX=192g      JAVA_XMS=4g
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
GRISMAN_CPUS=${GRISMAN_CPUS:-104}      # max out fennario
COMPSCI_CPUS=${COMPSCI_CPUS:-128}      # max out fitz
GRISMAN_MEM=${GRISMAN_MEM:-128G}
COMPSCI_MEM=${COMPSCI_MEM:-256G}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-4g}
TARGET=${TARGET:-both}                 # both | grisman | compsci

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
EXTRA_JVM_ARGS=${EXTRA_JVM_ARGS:-}
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx$JAVA_XMX -Xms$JAVA_XMS $EXTRA_JVM_ARGS"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

# Per-partition SLURM flags. grisman pins the 104-CPU a5000 (fennario) nodes;
# compsci asking for 128 cpus-per-task already restricts to the fitz nodes.
slurm_flags() {
    local part=$1 cpus=$2 mem=$3
    local f="--partition=$part --cpus-per-task=$cpus --mem=$mem --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL"
    if [ "$part" = "grisman" ]; then f="$f --account=grisman --constraint=a5000"; fi
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
        --wrap "$JAVA $JARGS \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=pac \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results \
            -Dosprey.bench.numCPUs=$cpus \
            -Dbranchmarkstar.usePAC=true \
            -Dbranchmarkstar.pac.samples=$PAC_SAMPLES \
            -Dbranchmarkstar.pac.confidence=$PAC_CONFIDENCE \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    echo "Submitted PAC $did ($pdb) [$part, ${cpus} CPUs, ${mem}, Xmx=${JAVA_XMX}]: $JID"
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
