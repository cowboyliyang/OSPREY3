#!/bin/bash
# One srun rank for sequence-bundle PACK* preflight.
# This script is intentionally not an sbatch script: the parent allocation
# launches one rank per node with srun.

set -euo pipefail

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
CP_FILE=${CP_FILE:-/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt}
: "${PDB_ROOT:?PDB_ROOT is required}"
: "${SRC:?SRC is required}"
: "${DESIGN_ID:?DESIGN_ID is required}"
: "${PDB:?PDB is required}"
: "${MUTABLE:?MUTABLE is required}"
: "${FLEXIBLE:?FLEXIBLE is required}"
: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${SHARDS:?SHARDS is required}"

# CP_FILE contains a few repository-relative entries (lib/*.jar).  Slurm
# preserves the submitter's working directory, so make their resolution
# independent of where sbatch was invoked.
if [ ! -d "$REPO" ]; then
    echo "REPO is not a directory: $REPO" >&2
    exit 2
fi
cd "$REPO"

RANK=${SLURM_PROCID:-0}
if [ "$SHARDS" -le 0 ] || [ "$RANK" -lt 0 ] || [ "$RANK" -ge "$SHARDS" ]; then
    echo "invalid rank=$RANK shards=$SHARDS" >&2
    exit 2
fi

THREADS=${SLURM_CPUS_PER_TASK:-16}
PREFLIGHT_ONLY=${PREFLIGHT_ONLY:-true}
GPU_COUNT_PER_NODE=${GPU_COUNT_PER_NODE:-${SLURM_GPUS_ON_NODE:-0}}
RANK_ROOT=$RESULT_ROOT/rank_$RANK
EDIR=$RANK_ROOT/emat_cache/$DESIGN_ID
POLICY_OUT=$RANK_ROOT/policy.tsv
SHARD_RESULT_OUT=$RANK_ROOT/shard_${RANK}.tsv
mkdir -p "$EDIR" "$RANK_ROOT/dp_mmap"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
fi
if [ "$PREFLIGHT_ONLY" != true ]; then
    if [ "$GPU_COUNT_PER_NODE" -le 0 ]; then
        echo "formal multi-node run requires visible GPUs; GPU_COUNT_PER_NODE=$GPU_COUNT_PER_NODE" >&2
        exit 3
    fi
    nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader || true
fi

# Each rank receives its own PACK* cache namespace. Symlinking the immutable
# export emats avoids making another large copy on the shared filesystem.
for st in complex ligand protein; do
    rigid="$EDIR/packstar.$st.rigid.dat"
    minimizing="$EDIR/packstar.$st.minimizing.dat"
    if [ ! -e "$rigid" ]; then
        ln -s "$SRC/export.$st.rigid.dat" "$rigid"
    fi
    if [ ! -e "$minimizing" ]; then
        ln -s "$SRC/export.$st.min.dat" "$minimizing"
    fi
done

if [ -n "${CLASS_ROOT:-}" ]; then
    # The checked-in classpath file contains the historical relative build
    # entries.  A Slurm validation build uses an xtmp buildDir, so rewrite only
    # those four entries while retaining the dependency jars from CP_FILE.
    CP=$(sed \
        -e "s#build/classes/java/main#${CLASS_ROOT}/classes/java/main#g" \
        -e "s#build/classes/java/test#${CLASS_ROOT}/classes/java/test#g" \
        -e "s#build/resources/main#${CLASS_ROOT}/resources/main#g" \
        -e "s#build/resources/test#${CLASS_ROOT}/resources/test#g" \
        "$CP_FILE")
else
    CP=$(<"$CP_FILE")
fi
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx80g -Xms4g -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$RANK_ROOT/dp_mmap"
    -Dpackstar.decomp.strategy=adaptive
    -Dpackstar.decomp.weightedHicks.restarts="${WEIGHTED_HICKS_RESTARTS:-32}"
    -Dpackstar.decomp.weightedHicks.randomMoves="${WEIGHTED_HICKS_RANDOM_MOVES:-32}"
    -Dpackstar.decomp.exactImprove.maxDrop="${EXACT_MAX_DROP:-2}"
    -Dpackstar.decomp.exactImprove.maxMillis="${INITIAL_EXACT_MILLIS:-120000}"
    -Dpackstar.pac.randomSeed="${PAC_RANDOM_SEED:-42}"
    -Dpackstar.rootSplit=predicted
    -Dpackstar.rootSplit.hostBudgetBytes="${HOST_BUDGET_BYTES:-450971566080}"
    -Dpackstar.rootSplit.gpuBudgetBytes="${GPU_BUDGET_BYTES:-23085449216}"
    -Dpackstar.dp.gpu=true
    -Dpackstar.dp.gpu.childSliceMaxBytes="${GPU_CHILD_SLICE_BYTES:-2147483648}"
    -Dpackstar.dp.gpu.outOfCore.outputWorkspaceMaxBytes="${GPU_OOC_WORKSPACE_BYTES:-4294967296}"
    -Dpackstar.admission.gpuWorkPerSecondPerGpu="${GPU_RATE_PER_GPU:-614000000}"
    -Dpackstar.admission.gpuCount="${GPU_COUNT:-8}"
    -Dpackstar.admission.oocBytesPerSecond="${OOC_RATE:-170000000}"
    -Dpackstar.admission.safetyFactor="${SAFETY_FACTOR:-1.5}"
    -Dpackstar.admission.softStateHours="${SOFT_STATE_HOURS:-24}"
    -Dpackstar.admission.caseSlaHours="${CASE_SLA_HOURS:-336}"
    -Dpackstar.admission.finalExactMaxMillis="${FINAL_EXACT_MILLIS:-180000}"
    -Dpackstar.admission.maxRounds="${MAX_ROUNDS:-4}"
    -Dpackstar.admission.previewThreads="$THREADS"
    -Dpackstar.admission.distributed=true
    -Dpackstar.admission.shardIndex="$RANK"
    -Dpackstar.admission.shardCount="$SHARDS"
    -Dpackstar.admission.shardSlaHours=1e300
    -Dpackstar.admission.forceFinalPass="${FORCE_FINAL_PASS:-true}"
    -Dpackstar.admission.policyOut="$POLICY_OUT"
    -Dpackstar.admission.shardResultOut="$SHARD_RESULT_OUT"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget="${BUDGET:-0.5}"
    -Dosprey.bench.method=packstar
    -Dosprey.bench.packstarPreflightOnly="$PREFLIGHT_ONLY"
    -Dosprey.bench.numCPUs="$THREADS"
    -Dosprey.bench.outputDir="$RANK_ROOT"
    -Dosprey.bench.designId="$DESIGN_ID"
    -Dosprey.bench.pdbPath="$PDB_ROOT/$PDB/$PDB.min.reduce.renum.pdb"
    "-Dosprey.bench.mutable=$MUTABLE"
    "-Dosprey.bench.flexible=$FLEXIBLE"
    -Dosprey.kstar.sequenceShardIndex="$RANK"
    -Dosprey.kstar.sequenceShardCount="$SHARDS"
)

if [ "$PREFLIGHT_ONLY" != true ]; then
    JARGS+=(
        -Dbranchdp.dp.gpu.multiGpu=true
        -Dbranchdp.dp.gpu.maxGpus="$GPU_COUNT_PER_NODE"
        -Dbranchdp.dp.gpu.minWork="${GPU_MIN_WORK:-0}"
        -Dbranchdp.dp.gpu.minMStatesPerGpu="${GPU_MIN_MSTATES_PER_GPU:-1}"
        -Dbranchdp.dp.gpu.failIfNoGpuPath=true
        -Dbranchdp.pac.sampling.gpu.multiGpu=true
        -Dbranchdp.pac.sampling.gpu.maxGpus="$GPU_COUNT_PER_NODE"
        -Dosprey.wmb.numGpus="$GPU_COUNT_PER_NODE"
        -Dosprey.wmb.streamsPerGpu="${GPU_STREAMS_PER_GPU:-8}"
    )
fi

if [ -n "${FINAL_MAX_STATES:-}" ]; then
    JARGS+=("-Dpackstar.admission.finalMaxStates=$FINAL_MAX_STATES")
fi

if [ -n "${POLICY_IN:-}" ]; then
    JARGS+=("-Dpackstar.admission.policyIn=$POLICY_IN")
fi
if [ "${LOCKED_POLICY_ONLY:-false}" = true ]; then
    if [ -z "${POLICY_IN:-}" ]; then
        echo "LOCKED_POLICY_ONLY=true requires POLICY_IN" >&2
        exit 2
    fi
    JARGS+=("-Dpackstar.admission.lockedPolicyOnly=true")
fi

MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
echo "pf_sequence_rank design=$DESIGN_ID pdb=$PDB budget=${BUDGET:-0.5} rank=$RANK/$SHARDS threads=$THREADS gpus=$GPU_COUNT_PER_NODE preflightOnly=$PREFLIGHT_ONLY node=${SLURMD_NODENAME:-unknown} start=$(date --iso-8601=seconds)"
echo "mutable=$MUTABLE"
echo "flexible=$FLEXIBLE"

"$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"

# The Java process writes the durable policy and shard TSV. Keep only those
# small audit products; the per-rank mmap directory is an intermediate.
rm -rf "$RANK_ROOT/dp_mmap"
echo "pf_sequence_rank_done design=$DESIGN_ID rank=$RANK time=$(date --iso-8601=seconds)"
