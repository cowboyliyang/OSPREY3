#!/bin/bash
# Formal PACK* K* single-node control run using the same worker and JVM
# settings as run_packstar_multigpu_multinode.slurm, but with one jerry node.
# This is a four-GPU single-node control for the two-node/four-GPU-per-node
# sequence-sharded validation.
#
#SBATCH --job-name=packstar_gpu_1n
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=jerry1
#SBATCH --nodes=1
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
[ -n "$design_id" ] || { echo "no design resolved for DESIGN=$DESIGN" >&2; exit 2; }

PDB_ROOT=${PDB_ROOT:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
SRC=${SRC:-/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id}
RESULT_ROOT=${RESULT_ROOT:-/usr/xtmp/lz280/packstar_production/formal_singlenode_gpu/${design_id}_rb${BUDGET}_${SLURM_JOB_ID}}
SHARDS=1
SHARD_COUNT=1
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-48}}

export REPO CLASS_ROOT RESULT_ROOT SHARDS SHARD_COUNT THREADS CSV DESIGN BUDGET
export DESIGN_ID="$design_id" PDB="$pdb" MUTABLE="$mutable" FLEXIBLE="$flexible"
export PDB_ROOT SRC
export PREFLIGHT_ONLY=false
export GPU_COUNT_PER_NODE=${GPU_COUNT_PER_NODE:-4}
export GPU_COUNT=${GPU_COUNT:-4}
export GPU_MIN_WORK=${GPU_MIN_WORK:-0}
export GPU_MIN_MSTATES_PER_GPU=${GPU_MIN_MSTATES_PER_GPU:-1}
export GPU_STREAMS_PER_GPU=${GPU_STREAMS_PER_GPU:-8}
export MAX_ROUNDS=${MAX_ROUNDS:-1}
export FINAL_MAX_STATES=${FINAL_MAX_STATES:-1}
export ADMISSION_METRIC=${ADMISSION_METRIC:-serial}
export CASE_SLA_HOURS=${CASE_SLA_HOURS:-336}
export SOFT_STATE_HOURS=${SOFT_STATE_HOURS:-24}

mkdir -p "$RESULT_ROOT"
echo "packstar GPU single-node control start job=$SLURM_JOB_ID node=$SLURM_NODELIST shards=1 gpusPerNode=$GPU_COUNT_PER_NODE threadsPerNode=$THREADS design=$DESIGN budget=$BUDGET start=$(date --iso-8601=seconds)"

set +e
srun --ntasks=1 --ntasks-per-node=1 --cpus-per-task="$THREADS" \
    --cpu-bind=cores bash "$REPO/slurm/scripts/packstar_preflight_sequence_rank.sh"
rank_status=$?
set -e

if [ "$rank_status" -ne 0 ]; then
    echo "single-node rank failure: srunStatus=$rank_status" >&2
    exit "$rank_status"
fi

export SLA_HOURS="$CASE_SLA_HOURS"
export REPORT_OUT="$RESULT_ROOT/summary.tsv"
export DETAIL_OUT="$RESULT_ROOT/detail.tsv"
export MERGED_POLICY_OUT="$RESULT_ROOT/merged_policy.tsv"
bash "$REPO/slurm/scripts/aggregate_packstar_preflight_sequence.sh"
export PACKSTAR_ALLOW_INCOMPLETE_RESULTS=${PACKSTAR_ALLOW_INCOMPLETE_RESULTS:-false}
bash "$REPO/slurm/scripts/merge_packstar_multinode_results.sh"

if [ "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS" = true ]; then
    echo "packstar GPU single-node diagnostic completed with explicit incomplete-result opt-in resultRoot=$RESULT_ROOT time=$(date --iso-8601=seconds)"
else
    echo "packstar GPU single-node production validation passed resultRoot=$RESULT_ROOT time=$(date --iso-8601=seconds)"
fi
