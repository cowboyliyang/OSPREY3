#!/bin/bash
# Sequence-bundle PACK* preflight for one (design, residualBudget) case.
# Submit as a plain job with DESIGN/BUDGET or as an array 0..113.
#
# The default admission metric remains serial because the current formal K*
# runner is single-node. Set ADMISSION_METRIC=makespan only when the matching
# sequence-bundle formal scheduler is being used. Both metrics are always
# written to the aggregate report.
#SBATCH --job-name=pf_seq_mn
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err

set -euo pipefail

CSV=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv
BUDGETS=(0.5 1 2)

if [ -n "${DESIGN:-}" ] && [ -n "${BUDGET:-}" ]; then
    row=$(awk -F, -v d="$DESIGN" 'NR>2 && $1==d {print; exit}' "$CSV")
    budget="$BUDGET"
else
    tid=${SLURM_ARRAY_TASK_ID:?need DESIGN+BUDGET or SLURM_ARRAY_TASK_ID}
    didx=$((tid/3))
    bidx=$((tid%3))
    budget=${BUDGETS[$bidx]}
    row=$(awk -F, -v n=$((didx+1)) 'NR>2 {c++; if(c==n){print; exit}}' "$CSV")
fi

design_id=$(printf '%s' "$row" | awk -F, '{print $1}')
pdb=$(printf '%s' "$row" | awk -F, '{print $2}')
mutable=$(printf '%s' "$row" | awk -F, '{print $6}')
flexible=$(printf '%s' "$row" | awk -F, '{print $7}')
[ -n "$design_id" ] || { echo "no design resolved" >&2; exit 2; }

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
CP_FILE=${CP_FILE:-/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt}
PDB_ROOT=${PDB_ROOT:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
SRC=${SRC:-/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id}
RESULT_ROOT=${RESULT_ROOT:-/usr/xtmp/lz280/packstar_production/preflight_multinode/${design_id}_rb${budget}_${SLURM_JOB_ID}}
SHARDS=${SHARDS:-${SLURM_NNODES:-1}}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
ADMISSION_METRIC=${ADMISSION_METRIC:-serial}

if [ "$SHARDS" -le 1 ]; then
    echo "multi-node preflight requires at least two nodes; SHARDS=$SHARDS" >&2
    exit 2
fi
if [ "$SHARDS" -gt "${SLURM_NNODES:-1}" ]; then
    echo "SHARDS=$SHARDS exceeds allocated nodes=${SLURM_NNODES:-1}" >&2
    exit 2
fi

for st in complex ligand protein; do
    [ -f "$SRC/export.$st.rigid.dat" ] || { echo "missing $SRC/export.$st.rigid.dat" >&2; exit 3; }
    [ -f "$SRC/export.$st.min.dat" ] || { echo "missing $SRC/export.$st.min.dat" >&2; exit 3; }
done

export REPO JAVA CP_FILE PDB_ROOT SRC DESIGN_ID="$design_id" PDB="$pdb"
export CLASS_ROOT=${CLASS_ROOT:-}
export MUTABLE="$mutable" FLEXIBLE="$flexible" RESULT_ROOT SHARDS THREADS
export SHARD_COUNT="$SHARDS"
export BUDGET="$budget"
export CASE_SLA_HOURS=${CASE_SLA_HOURS:-336}
export SOFT_STATE_HOURS=${SOFT_STATE_HOURS:-24}
export SAFETY_FACTOR=${SAFETY_FACTOR:-1.5}
export GPU_RATE_PER_GPU=${GPU_RATE_PER_GPU:-614000000}
export GPU_COUNT=${GPU_COUNT:-8}
export OOC_RATE=${OOC_RATE:-170000000}
export GPU_BUDGET_BYTES=${GPU_BUDGET_BYTES:-23085449216}
export HOST_BUDGET_BYTES=${HOST_BUDGET_BYTES:-450971566080}
export GPU_CHILD_SLICE_BYTES=${GPU_CHILD_SLICE_BYTES:-2147483648}
export GPU_OOC_WORKSPACE_BYTES=${GPU_OOC_WORKSPACE_BYTES:-4294967296}
export INITIAL_EXACT_MILLIS=${INITIAL_EXACT_MILLIS:-120000}
export FINAL_EXACT_MILLIS=${FINAL_EXACT_MILLIS:-180000}
export MAX_ROUNDS=${MAX_ROUNDS:-4}
export EXACT_MAX_DROP=${EXACT_MAX_DROP:-2}
export WEIGHTED_HICKS_RESTARTS=${WEIGHTED_HICKS_RESTARTS:-32}
export WEIGHTED_HICKS_RANDOM_MOVES=${WEIGHTED_HICKS_RANDOM_MOVES:-32}
export ADMISSION_METRIC

mkdir -p "$RESULT_ROOT"
echo "packstar sequence-bundle preflight start design=$design_id pdb=$pdb budget=$budget shards=$SHARDS threads=$THREADS metric=$ADMISSION_METRIC job=$SLURM_JOB_ID start=$(date --iso-8601=seconds)"
echo "mutable=$mutable"
echo "flexible=$flexible"
echo "resultRoot=$RESULT_ROOT"

set +e
srun --ntasks="$SHARDS" --ntasks-per-node=1 --cpus-per-task="$THREADS" \
    --cpu-bind=cores bash "$REPO/slurm/scripts/packstar_preflight_sequence_rank.sh"
rank_status=$?
set -e

if [ "$rank_status" -ne 0 ]; then
    echo "one or more sequence-bundle ranks failed, srunStatus=$rank_status" >&2
fi

export SLA_HOURS="$CASE_SLA_HOURS"
export REPORT_OUT="$RESULT_ROOT/summary.tsv"
export DETAIL_OUT="$RESULT_ROOT/detail.tsv"
export MERGED_POLICY_OUT="$RESULT_ROOT/merged_policy.tsv"

if [ "$rank_status" -eq 0 ]; then
    set +e
    bash "$REPO/slurm/scripts/aggregate_packstar_preflight_sequence.sh"
    aggregate_status=$?
    set -e
else
    aggregate_status=1
fi

echo "packstar sequence-bundle preflight end design=$design_id budget=$budget srunStatus=$rank_status aggregateStatus=$aggregate_status time=$(date --iso-8601=seconds)"
if [ "$rank_status" -ne 0 ]; then
    exit "$rank_status"
fi
exit "$aggregate_status"
