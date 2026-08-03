#!/bin/bash
# Production PACK* K* runner for one or more Titan V nodes.
#
# The allocation is used in three ordered phases:
#   1. a canonical one-shard preflight chooses a shard-invariant exact policy;
#   2. every requested shard previews that frozen policy, then the parent
#      aggregates a real global makespan and admits/rejects the case;
#   3. formal ranks load the admitted policy without repeating policy search,
#      followed by strict ownership, convergence, WT, and global-rank merging.
#
# Set FROZEN_POLICY_IN to reuse the exact same canonical policy across a
# single/2-shard/4-shard/repeat validation matrix.  Otherwise phase 1 creates
# a policy snapshot inside this job's result root.
#
#SBATCH --job-name=packstar_gpu_prod
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=grisman-40,jerry[1-2,4-7]
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --gres=gpu:titan_v:4
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
CLASS_ROOT=${CLASS_ROOT:-/usr/xtmp/lz280/osprey3-gradle-build/osprey-root}
CSV=${CSV:-/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv}
DESIGN=${DESIGN:-1a0r}
BUDGET=${BUDGET:-0.5}

row=$(awk -F, -v d="$DESIGN" 'NR>2 && $1==d {print; exit}' "$CSV")
design_id=$(printf '%s' "$row" | awk -F, '{print $1}')
pdb=$(printf '%s' "$row" | awk -F, '{print $2}')
mutable=$(printf '%s' "$row" | awk -F, '{print $6}')
flexible=$(printf '%s' "$row" | awk -F, '{print $7}')
[ -n "$design_id" ] || {
    echo "no design resolved for DESIGN=$DESIGN from CSV=$CSV" >&2
    exit 2
}

PDB_ROOT=${PDB_ROOT:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
SRC=${SRC:-/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id}
RESULT_ROOT=${RESULT_ROOT:-/usr/xtmp/lz280/packstar_production/formal_production_gpu/${design_id}_rb${BUDGET}_${SLURM_JOB_ID}}
SHARDS=${SHARDS:-${SLURM_NNODES:-1}}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-48}}
FROZEN_POLICY_IN=${FROZEN_POLICY_IN:-}
CASE_SLA_HOURS=${CASE_SLA_HOURS:-336}
POLICY_GENERATION_SLA_HOURS=${POLICY_GENERATION_SLA_HOURS:-1e300}
PACKSTAR_ALLOW_INCOMPLETE_RESULTS=${PACKSTAR_ALLOW_INCOMPLETE_RESULTS:-false}
CP_FILE_SOURCE=${CP_FILE:-/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt}

case "$RESULT_ROOT" in
    /usr/xtmp/lz280/*) ;;
    *)
        echo "RESULT_ROOT must be under /usr/xtmp/lz280, got $RESULT_ROOT" >&2
        exit 2
        ;;
esac
if ! [[ "$SHARDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "SHARDS must be a positive integer, got $SHARDS" >&2
    exit 2
fi
if [ -n "${SLURM_NNODES:-}" ] && [ "$SHARDS" -ne "$SLURM_NNODES" ]; then
    echo "SHARDS=$SHARDS must equal allocated SLURM_NNODES=$SLURM_NNODES" >&2
    exit 2
fi
case "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS" in
    true|false) ;;
    *)
        echo "PACKSTAR_ALLOW_INCOMPLETE_RESULTS must be true or false" >&2
        exit 2
        ;;
esac
if [ -n "$FROZEN_POLICY_IN" ] && [ ! -s "$FROZEN_POLICY_IN" ]; then
    echo "FROZEN_POLICY_IN is missing or empty: $FROZEN_POLICY_IN" >&2
    exit 2
fi
if [ ! -s "$CP_FILE_SOURCE" ]; then
    echo "classpath file is missing or empty: $CP_FILE_SOURCE" >&2
    exit 2
fi
if [ -e "$RESULT_ROOT/run_manifest.tsv" ] \
        || [ -d "$RESULT_ROOT/rank_0" ] \
        || [ -e "$RESULT_ROOT/frozen_policy.tsv" ]; then
    echo "refusing to overwrite an existing PACK* production result root: $RESULT_ROOT" >&2
    exit 2
fi

CANONICAL_ROOT=$RESULT_ROOT/canonical_preflight
DISTRIBUTED_ROOT=$RESULT_ROOT/distributed_preflight
FROZEN_POLICY=$RESULT_ROOT/frozen_policy.tsv
SCRIPT_SNAPSHOT_ROOT=$RESULT_ROOT/code_snapshot
mkdir -p "$SCRIPT_SNAPSHOT_ROOT/slurm/scripts" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/wrappers"
cp "$REPO/slurm/scripts/packstar_preflight_sequence_rank.sh" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/"
cp "$REPO/slurm/scripts/aggregate_packstar_preflight_sequence.sh" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/"
cp "$REPO/slurm/scripts/merge_packstar_multinode_results.sh" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/"
cp "$REPO/slurm/scripts/run_packstar_multinode_production.slurm" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/"
cp "$REPO/slurm/wrappers/packstar_multinode_merge.py" \
    "$SCRIPT_SNAPSHOT_ROOT/slurm/wrappers/"
cp "$CP_FILE_SOURCE" "$SCRIPT_SNAPSHOT_ROOT/classpath.txt"
WORKER_SCRIPT=$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/packstar_preflight_sequence_rank.sh
AGGREGATE_SCRIPT=$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/aggregate_packstar_preflight_sequence.sh
MERGE_SCRIPT=$SCRIPT_SNAPSHOT_ROOT/slurm/scripts/merge_packstar_multinode_results.sh
CP_FILE=$SCRIPT_SNAPSHOT_ROOT/classpath.txt

export REPO CLASS_ROOT THREADS CSV DESIGN BUDGET CP_FILE
export DESIGN_ID="$design_id" PDB="$pdb" MUTABLE="$mutable" FLEXIBLE="$flexible"
export PDB_ROOT SRC
export GPU_COUNT_PER_NODE=${GPU_COUNT_PER_NODE:-4}
export GPU_COUNT=${GPU_COUNT:-4}
export GPU_MIN_WORK=${GPU_MIN_WORK:-0}
export GPU_MIN_MSTATES_PER_GPU=${GPU_MIN_MSTATES_PER_GPU:-1}
export GPU_STREAMS_PER_GPU=${GPU_STREAMS_PER_GPU:-8}
export MAX_ROUNDS=${MAX_ROUNDS:-1}
export FINAL_MAX_STATES=${FINAL_MAX_STATES:-1}
export CASE_SLA_HOURS SOFT_STATE_HOURS=${SOFT_STATE_HOURS:-24}
export PAC_RANDOM_SEED=${PAC_RANDOM_SEED:-42}
export PACKSTAR_ALLOW_INCOMPLETE_RESULTS

run_rank_phase() (
    phase_name=$1
    phase_root=$2
    phase_shards=$3
    phase_preflight_only=$4
    phase_policy_in=$5
    phase_locked_policy_only=$6
    phase_force_final_pass=$7

    export RESULT_ROOT="$phase_root"
    export SHARDS="$phase_shards"
    export SHARD_COUNT="$phase_shards"
    export PREFLIGHT_ONLY="$phase_preflight_only"
    export POLICY_IN="$phase_policy_in"
    export LOCKED_POLICY_ONLY="$phase_locked_policy_only"
    export FORCE_FINAL_PASS="$phase_force_final_pass"

    echo "packstar production phase=$phase_name shards=$phase_shards preflightOnly=$phase_preflight_only lockedPolicyOnly=$phase_locked_policy_only forceFinalPass=$phase_force_final_pass root=$phase_root start=$(date --iso-8601=seconds)"
    srun --nodes="$phase_shards" --ntasks="$phase_shards" \
        --ntasks-per-node=1 --cpus-per-task="$THREADS" \
        --cpu-bind=cores --kill-on-bad-exit=1 \
        bash "$WORKER_SCRIPT"
    echo "packstar production phase=$phase_name complete time=$(date --iso-8601=seconds)"
)

aggregate_phase() (
    phase_name=$1
    phase_root=$2
    phase_shards=$3
    phase_metric=$4
    phase_sla_hours=$5

    export RESULT_ROOT="$phase_root"
    export SHARD_COUNT="$phase_shards"
    export ADMISSION_METRIC="$phase_metric"
    export SLA_HOURS="$phase_sla_hours"
    export REPORT_OUT="$phase_root/summary.tsv"
    export DETAIL_OUT="$phase_root/detail.tsv"
    export MERGED_POLICY_OUT="$phase_root/merged_policy.tsv"

    echo "packstar production aggregate=$phase_name shards=$phase_shards metric=$phase_metric slaHours=$phase_sla_hours"
    bash "$AGGREGATE_SCRIPT"
)

echo "packstar production start job=$SLURM_JOB_ID nodes=$SLURM_NODELIST shards=$SHARDS gpusPerNode=$GPU_COUNT_PER_NODE threadsPerNode=$THREADS design=$design_id budget=$BUDGET seed=$PAC_RANDOM_SEED start=$(date --iso-8601=seconds)"

if [ -z "$FROZEN_POLICY_IN" ]; then
    # Policy selection always sees the whole sequence set in one canonical
    # shard.  Its deliberately unbounded aggregation SLA is only for choosing
    # the plan; the real allocation is admitted in the next phase.
    run_rank_phase canonical-policy "$CANONICAL_ROOT" 1 true "" false true
    aggregate_phase canonical-policy "$CANONICAL_ROOT" 1 serial \
        "$POLICY_GENERATION_SLA_HOURS"
    cp "$CANONICAL_ROOT/merged_policy.tsv" "$FROZEN_POLICY"
    policy_source=generated
else
    cp "$FROZEN_POLICY_IN" "$FROZEN_POLICY"
    policy_source=$FROZEN_POLICY_IN
fi

if [ ! -s "$FROZEN_POLICY" ]; then
    echo "frozen policy snapshot is empty: $FROZEN_POLICY" >&2
    exit 3
fi
policy_sha256=$(sha256sum "$FROZEN_POLICY" | awk '{print $1}')

# Every actual shard now previews only the one frozen policy.  forceFinalPass
# must be false in lock mode: a forced final pass means policy search, which is
# intentionally forbidden after the canonical policy has been selected.
run_rank_phase distributed-admission "$DISTRIBUTED_ROOT" "$SHARDS" true \
    "$FROZEN_POLICY" false false
aggregate_phase distributed-admission "$DISTRIBUTED_ROOT" "$SHARDS" \
    makespan "$CASE_SLA_HOURS"

{
    printf '%b\n' 'schema\tpackstar-production-run-v1'
    printf 'jobId\t%s\n' "$SLURM_JOB_ID"
    printf 'designId\t%s\n' "$design_id"
    printf 'residualBudget\t%s\n' "$BUDGET"
    printf 'shardCount\t%s\n' "$SHARDS"
    printf 'nodes\t%s\n' "$SLURM_NODELIST"
    printf 'threadsPerRank\t%s\n' "$THREADS"
    printf 'gpusPerRank\t%s\n' "$GPU_COUNT_PER_NODE"
    printf 'pacRandomSeed\t%s\n' "$PAC_RANDOM_SEED"
    printf 'classRoot\t%s\n' "$CLASS_ROOT"
    printf 'caseSlaHours\t%s\n' "$CASE_SLA_HOURS"
    printf 'policySource\t%s\n' "$policy_source"
    printf 'frozenPolicy\t%s\n' "$FROZEN_POLICY"
    printf 'frozenPolicySha256\t%s\n' "$policy_sha256"
    printf 'distributedPreflight\t%s\n' "$DISTRIBUTED_ROOT/detail.tsv"
    printf 'allowIncompleteResults\t%s\n' "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS"
} > "$RESULT_ROOT/run_manifest.tsv"

# The parent has already admitted the exact distributed manifest.  Formal
# ranks load the same policy and skip only the redundant policy-search preview;
# Branch-DP still reconstructs and checks each retained prediction ceiling
# before materializing its DP table.
run_rank_phase formal "$RESULT_ROOT" "$SHARDS" false \
    "$FROZEN_POLICY" true false

(
    export REPO="$SCRIPT_SNAPSHOT_ROOT"
    export RESULT_ROOT SHARD_COUNT="$SHARDS" DESIGN_ID="$design_id"
    export PREFLIGHT_DETAIL="$DISTRIBUTED_ROOT/detail.tsv"
    export PACKSTAR_ALLOW_INCOMPLETE_RESULTS
    bash "$MERGE_SCRIPT"
)

if [ "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS" = true ]; then
    echo "packstar production diagnostic complete (incomplete-result opt-in) resultRoot=$RESULT_ROOT time=$(date --iso-8601=seconds)"
else
    echo "packstar production validation passed resultRoot=$RESULT_ROOT time=$(date --iso-8601=seconds)"
fi
