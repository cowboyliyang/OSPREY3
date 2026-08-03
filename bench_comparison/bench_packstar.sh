#!/bin/bash
# ===================================================================
# PACK* PAC-mode benchmark on dance_bench PDBs
#
# PAC = Probably Approximately Correct partition-function estimation via
# Rao-Blackwellized importance sampling on the branch-decomposition DP tree.
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
# Usage: bash bench_packstar.sh [design_id or "all"]
#   env overrides: PACKSTAR_SAMPLES=1000  PACKSTAR_CONFIDENCE=0.05  PACKSTAR_RESIDUAL_BOUND=1.0
#                  PACKSTAR_TRAIN_SAMPLES= PACKSTAR_PILOT_SAMPLES= PACKSTAR_MAX_EST_SAMPLES=
#                  PACKSTAR_SAMPLING_BATCHED=true PACKSTAR_SAMPLING_GPU=false PACKSTAR_SAMPLING_THREADS= PACKSTAR_SAMPLING_LARGE_LAMBDA=65536
#                  PACKSTAR_CLIP=true PACKSTAR_CLIP_QUANTILE=0.85 PACKSTAR_ITERATE=true PACKSTAR_ITERATE_MAX_ROUNDS=4 PACKSTAR_SIZE_SAFETY=0.9
#                  ROOT_SPLIT=memory
#                  GRISMAN_CPUS=104  COMPSCI_CPUS=128
#                  GRISMAN_MEM=400G   COMPSCI_MEM=256G
#                  JAVA_XMX=192g      JAVA_XMS=4g
#                  DP_CACHE=false     EXTRA_JVM_ARGS="..."
#                  GRISMAN_FULL_NODE=false GRISMAN_EXCLUDE="" GRISMAN_CONSTRAINT=a5000
#                  (set GRISMAN_CONSTRAINT= to allow any grisman node)
#                  TARGET=both|grisman|compsci   (default: grisman)
# Compatibility: bench_pac.sh and old PAC_* env vars remain accepted for old commands.
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
DESIGN=${1:-2q1e}

# --- tunables ---
PACKSTAR_SAMPLES=${PACKSTAR_SAMPLES:-${PAC_SAMPLES:-1000}}       # packstar.pac.samples (budget: train 50% / pilot 10% / est rest)
PACKSTAR_CONFIDENCE=${PACKSTAR_CONFIDENCE:-${PAC_CONFIDENCE:-0.05}} # delta; 0.05 => 95% confidence
PACKSTAR_RESIDUAL_BOUND=${PACKSTAR_RESIDUAL_BOUND:-${PAC_RESIDUAL_BOUND:-1.0}} # deterministic |xi| bound, kcal/mol
PACKSTAR_TARGET_EPSILON=${PACKSTAR_TARGET_EPSILON:-${PAC_TARGET_EPSILON:-}}
PACKSTAR_TRAIN_SAMPLES=${PACKSTAR_TRAIN_SAMPLES:-${PAC_TRAIN_SAMPLES:-}}
PACKSTAR_PILOT_SAMPLES=${PACKSTAR_PILOT_SAMPLES:-${PAC_PILOT_SAMPLES:-}}
PACKSTAR_MAX_EST_SAMPLES=${PACKSTAR_MAX_EST_SAMPLES:-${PAC_MAX_EST_SAMPLES:-4000}}   # generous cap so hard seqs never undersize (6/7: N* up to 2305)
PACKSTAR_SAMPLING_BATCHED=${PACKSTAR_SAMPLING_BATCHED:-${PAC_SAMPLING_BATCHED:-true}}
PACKSTAR_SAMPLING_GPU=${PACKSTAR_SAMPLING_GPU:-${PAC_SAMPLING_GPU:-false}}
PACKSTAR_SAMPLING_THREADS=${PACKSTAR_SAMPLING_THREADS:-${PAC_SAMPLING_THREADS:-}}
PACKSTAR_SAMPLING_LARGE_LAMBDA=${PACKSTAR_SAMPLING_LARGE_LAMBDA:-${PAC_SAMPLING_LARGE_LAMBDA:-65536}}
PACKSTAR_CLIP=${PACKSTAR_CLIP:-${PAC_CLIP:-true}}              # packstar.pac.clip (6/7 behavior = true)
PACKSTAR_CLIP_QUANTILE=${PACKSTAR_CLIP_QUANTILE:-${PAC_CLIP_QUANTILE:-0.85}}
PACKSTAR_ITERATE=${PACKSTAR_ITERATE:-${PAC_ITERATE:-true}}
PACKSTAR_ITERATE_MAX_ROUNDS=${PACKSTAR_ITERATE_MAX_ROUNDS:-${PAC_ITERATE_MAX_ROUNDS:-4}}
PACKSTAR_SIZE_SAFETY=${PACKSTAR_SIZE_SAFETY:-${PAC_SIZE_SAFETY:-0.9}}
GRISMAN_CPUS=${GRISMAN_CPUS:-104}      # max out fennario
COMPSCI_CPUS=${COMPSCI_CPUS:-128}      # max out fitz
GRISMAN_MEM=${GRISMAN_MEM:-400G}
COMPSCI_MEM=${COMPSCI_MEM:-256G}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-192g}   # commit heap (==Xmx) so big DP tables don't trigger GC-thrash that serializes the 104-thread CCD
DP_CACHE=${DP_CACHE:-false}
ROOT_SPLIT=${ROOT_SPLIT:-memory}
GRISMAN_FULL_NODE=${GRISMAN_FULL_NODE:-false}
GRISMAN_EXCLUDE=${GRISMAN_EXCLUDE:-}
GRISMAN_CONSTRAINT=${GRISMAN_CONSTRAINT-a5000}
DP_GPU=${DP_GPU:-true}                  # true => GPU full-DP fast path (dp.gpu)
GRISMAN_GPUS=${GRISMAN_GPUS:-8}         # non-empty => --gres=gpu:a5000:N on grisman
TARGET=${TARGET:-grisman}               # both | grisman | compsci
EXCLUDE=${EXCLUDE:-3bua 3k3q 4z80}      # space-separated design ids to skip (for "all")

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord -XX:+UseParallelGC -Xlog:gc::uptime,level,tags"
if [ "$DP_GPU" = "true" ]; then
    EXTRA_JVM_ARGS="$EXTRA_JVM_ARGS -Dpackstar.dp.gpu=true -Dpackstar.dp.gpu.multiGpu=true -Dpackstar.dp.gpu.persistentContext=true"
fi
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
        if [ -n "$GRISMAN_GPUS" ]; then f="$f --gres=gpu:a5000:$GRISMAN_GPUS"; fi
    else
        f="$f --cpus-per-task=$cpus --mem=$mem"
    fi
    echo "$f"
}

mkdir -p "$LOGDIR" "$OUTDIR/results"

./gradlew testClasses --no-daemon -Dorg.gradle.vfs.watch=false 2>&1 | tail -3
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

    JID=$(sbatch $(slurm_flags "$part" "$cpus" "$mem") --job-name=packstar_${did} \
        --output=$LOGDIR/packstar_${did}_%j.out --error=$LOGDIR/packstar_${did}_%j.err \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-$cpus}; $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=packstar \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results \
            -Dosprey.bench.numCPUs=\$RUN_CPUS \
            -Dpackstar.pac.residualBound=$PACKSTAR_RESIDUAL_BOUND \
            -Dpackstar.rootSplit=$ROOT_SPLIT \
            -Dpackstar.dp.cache=$DP_CACHE \
            -Dpackstar.pac.samples=$PACKSTAR_SAMPLES \
            -Dpackstar.pac.confidence=$PACKSTAR_CONFIDENCE \
            -Dpackstar.pac.targetEpsilon=$PACKSTAR_TARGET_EPSILON \
            -Dpackstar.pac.trainSamples=$PACKSTAR_TRAIN_SAMPLES \
            -Dpackstar.pac.pilotSamples=$PACKSTAR_PILOT_SAMPLES \
            -Dpackstar.pac.maxEstSamples=$PACKSTAR_MAX_EST_SAMPLES \
            -Dpackstar.pac.sampling.batched=$PACKSTAR_SAMPLING_BATCHED \
            -Dpackstar.pac.sampling.gpu=$PACKSTAR_SAMPLING_GPU \
            -Dpackstar.pac.sampling.threads=$PACKSTAR_SAMPLING_THREADS \
            -Dpackstar.pac.sampling.largeLambdaThreshold=$PACKSTAR_SAMPLING_LARGE_LAMBDA \
            -Dpackstar.pac.clip=$PACKSTAR_CLIP \
            -Dpackstar.pac.clipQuantile=$PACKSTAR_CLIP_QUANTILE \
            -Dpackstar.pac.iterate=$PACKSTAR_ITERATE \
            -Dpackstar.pac.iterate.maxRounds=$PACKSTAR_ITERATE_MAX_ROUNDS \
            -Dpackstar.pac.sizeSafety=$PACKSTAR_SIZE_SAFETY \
            -cp \"\$(cat $LOGDIR/.classpath_bench.txt)\" $MAIN 2>&1" \
        | awk '{print $4}')
    if [ "$part" = "grisman" ] && [ "$GRISMAN_FULL_NODE" = "true" ]; then
        echo "Submitted PACK* $did ($pdb) [$part, full-node, exclude=${GRISMAN_EXCLUDE:-none}, constraint=${GRISMAN_CONSTRAINT:-none}, Xmx=${JAVA_XMX}]: $JID"
    else
        echo "Submitted PACK* $did ($pdb) [$part, ${cpus} CPUs, ${mem}, Xmx=${JAVA_XMX}]: $JID"
    fi
}

if [ "$DESIGN" = "all" ]; then
    i=0
    for did in $(grep -v "^#" "$SPECS" | grep -v "^[[:space:]]*$" | cut -d',' -f1); do
        skip=false
        for ex in $EXCLUDE; do if [ "$did" = "$ex" ]; then skip=true; break; fi; done
        if [ "$skip" = "true" ]; then echo "Skipping excluded design: $did"; continue; fi
        case "$TARGET" in
            grisman) submit_one "$did" grisman "$GRISMAN_CPUS" ;;
            compsci) submit_one "$did" compsci "$COMPSCI_CPUS" ;;
            both)    if [ $((i % 2)) -eq 0 ]; then submit_one "$did" grisman "$GRISMAN_CPUS";
                     else submit_one "$did" compsci "$COMPSCI_CPUS"; fi ;;
        esac
        i=$((i + 1))
    done
else
    # one or more designs (space-separated): default to grisman unless TARGET=compsci
    for did in $DESIGN; do
        if [ "$TARGET" = "compsci" ]; then submit_one "$did" compsci "$COMPSCI_CPUS";
        else submit_one "$did" grisman "$GRISMAN_CPUS"; fi
    done
fi
