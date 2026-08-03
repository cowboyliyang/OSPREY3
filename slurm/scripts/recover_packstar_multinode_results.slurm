#!/bin/bash
# Re-run only preflight aggregation and canonical result validation for a
# completed set of rank artifacts.  No PACK* estimator or DP computation runs.

#SBATCH --job-name=packstar_recover
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
export REPO
: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${SHARD_COUNT:?SHARD_COUNT is required}"
: "${DESIGN_ID:?DESIGN_ID is required}"

export SLA_HOURS=${SLA_HOURS:-336}
export ADMISSION_METRIC=${ADMISSION_METRIC:-makespan}
export REPORT_OUT=${REPORT_OUT:-$RESULT_ROOT/summary.tsv}
export DETAIL_OUT=${DETAIL_OUT:-$RESULT_ROOT/detail.tsv}
export MERGED_POLICY_OUT=${MERGED_POLICY_OUT:-$RESULT_ROOT/merged_policy.tsv}
export PACKSTAR_ALLOW_INCOMPLETE_RESULTS=${PACKSTAR_ALLOW_INCOMPLETE_RESULTS:-false}

bash "$REPO/slurm/scripts/aggregate_packstar_preflight_sequence.sh"
bash "$REPO/slurm/scripts/merge_packstar_multinode_results.sh"

echo "packstar result recovery complete resultRoot=$RESULT_ROOT time=$(date --iso-8601=seconds)"
