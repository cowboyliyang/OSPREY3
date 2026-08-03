#!/bin/bash
# Validate rank CSV ownership/convergence and publish one canonical global CSV.
# Run only inside the parent Slurm allocation after preflight aggregation.

set -euo pipefail

: "${REPO:?REPO is required}"
: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${SHARD_COUNT:?SHARD_COUNT is required}"
: "${DESIGN_ID:?DESIGN_ID is required}"

PYTHON=${PYTHON:-python3}
PACKSTAR_ALLOW_INCOMPLETE_RESULTS=${PACKSTAR_ALLOW_INCOMPLETE_RESULTS:-false}
MERGED_RESULT_OUT=${MERGED_RESULT_OUT:-$RESULT_ROOT/merged_packstar.csv}
INCOMPLETE_RESULT_OUT=${INCOMPLETE_RESULT_OUT:-$RESULT_ROOT/merged_packstar.incomplete.csv}
RESULT_VALIDATION_REPORT=${RESULT_VALIDATION_REPORT:-$RESULT_ROOT/result_validation.json}
RESULT_PROVENANCE_OUT=${RESULT_PROVENANCE_OUT:-$RESULT_ROOT/result_provenance.tsv}
PREFLIGHT_DETAIL=${PREFLIGHT_DETAIL:-$RESULT_ROOT/detail.tsv}

case "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS" in
    true|false) ;;
    *)
        echo "PACKSTAR_ALLOW_INCOMPLETE_RESULTS must be true or false" >&2
        exit 2
        ;;
esac

args=(
    --result-root "$RESULT_ROOT"
    --detail "$PREFLIGHT_DETAIL"
    --shard-count "$SHARD_COUNT"
    --design-id "$DESIGN_ID"
    --output "$MERGED_RESULT_OUT"
    --incomplete-output "$INCOMPLETE_RESULT_OUT"
    --report "$RESULT_VALIDATION_REPORT"
    --provenance "$RESULT_PROVENANCE_OUT"
)
if [ "$PACKSTAR_ALLOW_INCOMPLETE_RESULTS" = true ]; then
    args+=(--allow-incomplete)
fi

"$PYTHON" "$REPO/slurm/wrappers/packstar_multinode_merge.py" "${args[@]}"
