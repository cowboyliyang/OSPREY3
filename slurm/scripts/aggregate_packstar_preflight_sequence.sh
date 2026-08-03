#!/bin/bash
# Reduce sequence-bundle preflight products from one multi-node allocation.
# Run this inside the parent Slurm allocation, after all srun ranks finish.

set -euo pipefail

: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${SHARD_COUNT:?SHARD_COUNT is required}"

SLA_HOURS=${SLA_HOURS:-336}
ADMISSION_METRIC=${ADMISSION_METRIC:-serial}
REPORT_OUT=${REPORT_OUT:-$RESULT_ROOT/summary.tsv}
DETAIL_OUT=${DETAIL_OUT:-$RESULT_ROOT/detail.tsv}
MERGED_POLICY_OUT=${MERGED_POLICY_OUT:-$RESULT_ROOT/merged_policy.tsv}

case "$ADMISSION_METRIC" in
    serial|makespan) ;;
    *) echo "ADMISSION_METRIC must be serial or makespan, got $ADMISSION_METRIC" >&2; exit 2 ;;
esac

if [ "$SHARD_COUNT" -le 0 ]; then
    echo "SHARD_COUNT must be positive" >&2
    exit 2
fi

mkdir -p "$RESULT_ROOT"
DETAIL_TMP=$DETAIL_OUT.tmp.$$
SUMMARY_TMP=$REPORT_OUT.tmp.$$
POLICY_ROWS_TMP=$MERGED_POLICY_OUT.rows.tmp.$$
POLICY_TMP=$MERGED_POLICY_OUT.tmp.$$
cleanup() {
    rm -f "$DETAIL_TMP" "$SUMMARY_TMP" "$POLICY_ROWS_TMP" "$POLICY_TMP"
}
trap cleanup EXIT

printf '%b\n' 'shardIndex\tbundleOrdinal\tglobalSequenceB64\tordinal\tstateNameB64\tsequenceB64\tstateKeyB64\tbranchwidth\trootSplitEdge\tpredictedSeconds\tgpuWork\toocTrafficBytes\toocTrafficAvailable\tdpSweeps\tadaptiveAttempted\tadaptiveAccepted' > "$DETAIL_TMP"

TOTAL_BUNDLES=
TOTAL_STATES=
REPLICATED_BUNDLE_ORDINAL=
TOTAL_ROWS=0
TOTAL_SCHEDULED_BUNDLES=0
SERIAL_SECONDS=0
MAX_SHARD_SECONDS=0

meta_value() {
    local key=$1
    local file=$2
    sed -n "s/^# ${key}=//p" "$file" | head -n 1
}

for ((rank=0; rank<SHARD_COUNT; rank++)); do
    file=$RESULT_ROOT/rank_${rank}/shard_${rank}.tsv
    if [ ! -s "$file" ]; then
        echo "missing shard result: $file" >&2
        exit 3
    fi

    shard_index=$(meta_value shardIndex "$file")
    shard_count=$(meta_value shardCount "$file")
    total_bundles=$(meta_value totalBundles "$file")
    total_states=$(meta_value totalStates "$file")
    states_per_bundle=$(meta_value statesPerBundle "$file")
    assigned_bundles=$(meta_value assignedBundles "$file")
    replicated_bundle_ordinal=$(meta_value replicatedBundleOrdinal "$file")
    if [ "$shard_index" != "$rank" ] || [ "$shard_count" != "$SHARD_COUNT" ]; then
        echo "shard metadata mismatch in $file: index=$shard_index count=$shard_count" >&2
        exit 3
    fi
    if ! [[ "$total_bundles" =~ ^[0-9]+$ && "$total_states" =~ ^[0-9]+$ \
        && "$states_per_bundle" =~ ^[0-9]+$ && "$assigned_bundles" =~ ^[0-9]+$ \
        && "$replicated_bundle_ordinal" =~ ^-?[0-9]+$ ]]; then
        echo "invalid bundle/state/replication metadata in $file" >&2
        exit 3
    fi
    if [ "$total_bundles" -gt 0 ] && [ "$states_per_bundle" -le 0 ]; then
        echo "statesPerBundle must be positive when bundles are present in $file" >&2
        exit 3
    fi
    if [ "$total_bundles" -gt 0 ] && [ "$total_states" -ne $((total_bundles * states_per_bundle)) ]; then
        echo "totalStates is inconsistent with totalBundles/statesPerBundle in $file" >&2
        exit 3
    fi
    if [ -z "$TOTAL_BUNDLES" ]; then
        TOTAL_BUNDLES=$total_bundles
        TOTAL_STATES=$total_states
    elif [ "$TOTAL_BUNDLES" != "$total_bundles" ] || [ "$TOTAL_STATES" != "$total_states" ]; then
        echo "global metadata mismatch in $file" >&2
        exit 3
    fi
    if [ -z "${STATES_PER_BUNDLE:-}" ]; then
        STATES_PER_BUNDLE=$states_per_bundle
    elif [ "$STATES_PER_BUNDLE" != "$states_per_bundle" ]; then
        echo "statesPerBundle mismatch in $file" >&2
        exit 3
    fi
    if [ -z "$REPLICATED_BUNDLE_ORDINAL" ]; then
        REPLICATED_BUNDLE_ORDINAL=$replicated_bundle_ordinal
    elif [ "$REPLICATED_BUNDLE_ORDINAL" != "$replicated_bundle_ordinal" ]; then
        echo "replicatedBundleOrdinal mismatch in $file" >&2
        exit 3
    fi
    if [ "$replicated_bundle_ordinal" -lt -1 ] \
        || [ "$replicated_bundle_ordinal" -ge "$total_bundles" ]; then
        echo "replicatedBundleOrdinal out of range in $file" >&2
        exit 3
    fi

    stats=$(awk -F '\t' -v rank="$rank" -v detail="$DETAIL_TMP" '
        BEGIN { rows=0; bundles=0; sum=0; seenHeader=0 }
        /^#/ { next }
        /^bundleOrdinal\t/ { seenHeader=1; next }
        {
            if (NF != 15) {
                printf("bad field count=%d at input line %d\n", NF, NR) > "/dev/stderr";
                exit 10;
            }
            if ($1 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {
                print "non-integer bundle/ordinal" > "/dev/stderr";
                exit 11;
            }
            if ($9 !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ || ($9 + 0) < 0) {
                print "non-finite or negative predictedSeconds" > "/dev/stderr";
                exit 12;
            }
            print rank "\t" $0 >> detail;
            rows++;
            sum += $9;
            if (!seen[$1]++) bundles++;
        }
        END {
            if (!seenHeader) {
                print "missing v3 shard header" > "/dev/stderr";
                exit 13;
            }
            printf("%.17g\t%d\t%d\n", sum, rows, bundles);
        }
    ' "$file")
    read -r shard_seconds shard_rows shard_bundles <<< "$stats"
    if [ "$shard_bundles" -ne "$assigned_bundles" ] \
        || [ "$shard_rows" -ne $((assigned_bundles * states_per_bundle)) ]; then
        echo "assigned bundle/state count mismatch in $file" >&2
        exit 3
    fi
    TOTAL_ROWS=$((TOTAL_ROWS + shard_rows))
    TOTAL_SCHEDULED_BUNDLES=$((TOTAL_SCHEDULED_BUNDLES + shard_bundles))
    SERIAL_SECONDS=$(awk -v a="$SERIAL_SECONDS" -v b="$shard_seconds" 'BEGIN { printf "%.17g", a+b }')
    MAX_SHARD_SECONDS=$(awk -v a="$MAX_SHARD_SECONDS" -v b="$shard_seconds" 'BEGIN { print (b>a ? b : a) }')
done

GLOBAL_STATS=$(awk -F '\t' -v totalBundles="$TOTAL_BUNDLES" -v totalStates="$TOTAL_STATES" \
    -v statesPerBundle="$STATES_PER_BUNDLE" -v shardCount="$SHARD_COUNT" \
    -v replicatedBundle="$REPLICATED_BUNDLE_ORDINAL" '
    NR == 1 { next }
    {
        shard=$1; bundle=$2; ordinal=$4; seconds=$10 + 0;
        if (shard !~ /^[0-9]+$/ || shard < 0 || shard >= shardCount) exit 20;
        if (bundle !~ /^[0-9]+$/ || ordinal !~ /^[0-9]+$/) exit 20;
        if (bundle < 0 || bundle >= totalBundles) exit 21;
        if (ordinal < 0 || ordinal >= totalStates) exit 22;
        if (seenBundle[bundle]++ == 0) bundleCount++;
        bundleRows[bundle]++;
        bundleShard[bundle SUBSEP shard]=1;
        if (seenShardOrdinal[shard SUBSEP ordinal]++) exit 23;
        ordinalRows[ordinal]++;
        payload=$2;
        for (i=3; i<=NF; i++) payload=payload "\t" $i;
        if (bundle == replicatedBundle) {
            if (replicatedPayload[ordinal] != "" && replicatedPayload[ordinal] != payload) exit 24;
            replicatedPayload[ordinal]=payload;
        } else if (ordinalRows[ordinal] > 1) {
            exit 25;
        }
        serial += seconds;
        perShard[shard] += seconds;
        rows++;
    }
    END {
        expectedRows=totalStates;
        if (replicatedBundle >= 0) {
            expectedRows += (shardCount - 1) * statesPerBundle;
        }
        if (rows != expectedRows || bundleCount != totalBundles) exit 26;
        for (b=0; b<totalBundles; b++) {
            expectedCopies=(b == replicatedBundle ? shardCount : 1);
            if (bundleRows[b] != statesPerBundle * expectedCopies) exit 27;
            shardCopies=0;
            for (s=0; s<shardCount; s++) {
                if (bundleShard[b SUBSEP s]) shardCopies++;
            }
            if (shardCopies != expectedCopies) exit 28;
        }
        for (o=0; o<totalStates; o++) {
            isReplicatedOrdinal=(replicatedBundle >= 0 && o >= replicatedBundle * statesPerBundle && o < (replicatedBundle + 1) * statesPerBundle);
            expectedOrdinalCopies=(isReplicatedOrdinal ? shardCount : 1);
            if (ordinalRows[o] != expectedOrdinalCopies) exit 29;
        }
        max=0;
        for (s in perShard) if (perShard[s] > max) max=perShard[s];
        printf("%.17g\t%.17g\t%d\t%d\n", serial, max, rows, bundleCount);
    }
' "$DETAIL_TMP") || {
    echo "global shard validation failed (mutants must have one owner; replicated WT must be identical on every shard)" >&2
    exit 4
}

read -r SERIAL_SECONDS CHECKED_MAX_SECONDS CHECKED_ROWS CHECKED_BUNDLES <<< "$GLOBAL_STATS"
if [ "$CHECKED_ROWS" -ne "$TOTAL_ROWS" ] || [ "$CHECKED_BUNDLES" -ne "$TOTAL_BUNDLES" ]; then
    echo "global row/bundle count mismatch" >&2
    exit 4
fi

# The global recomputation is authoritative; the per-file accumulation above is
# retained only as an early diagnostic for a missing or malformed rank.
MAX_SHARD_SECONDS=$CHECKED_MAX_SECONDS
SERIAL_HOURS=$(awk -v s="$SERIAL_SECONDS" 'BEGIN { printf "%.10f", s/3600.0 }')
MAKESPAN_HOURS=$(awk -v s="$MAX_SHARD_SECONDS" 'BEGIN { printf "%.10f", s/3600.0 }')
if [ "$ADMISSION_METRIC" = makespan ]; then
    ADMISSION_HOURS=$MAKESPAN_HOURS
else
    ADMISSION_HOURS=$SERIAL_HOURS
fi
ADMITTED=$(awk -v h="$ADMISSION_HOURS" -v sla="$SLA_HOURS" 'BEGIN { print (h <= sla ? "true" : "false") }')

# Merge exact policies while rejecting conflicting plans for the same state key.
# Identical duplicate rows are harmless because independent sequence bundles can
# contain the same filtered state on different ranks.
: > "$POLICY_ROWS_TMP"
for ((rank=0; rank<SHARD_COUNT; rank++)); do
    policy=$RESULT_ROOT/rank_${rank}/policy.tsv
    if [ ! -f "$policy" ]; then
        echo "missing policy file: $policy" >&2
        exit 5
    fi
    awk '/^#/ || NF==0 { next } { print }' "$policy" >> "$POLICY_ROWS_TMP"
done

awk -F '\t' '
    BEGIN { print "# packstar exact-policy merged dump v1: minDrop\tmaxDrop\tmaxMillis\tmaxPredictedSeconds\tstateKey" }
    {
        if (NF < 5) exit 30;
        key=$NF;
        if (key in seen) {
            if (seen[key] != $0) exit 31;
            next;
        }
        seen[key]=$0;
        print;
    }
' "$POLICY_ROWS_TMP" > "$POLICY_TMP" || {
    echo "conflicting or malformed exact policies" >&2
    exit 5
}
mv "$POLICY_TMP" "$MERGED_POLICY_OUT"
POLICY_ROWS=$(awk 'NF && $0 !~ /^#/ { n++ } END { print n+0 }' "$MERGED_POLICY_OUT")

{
    printf '%b\n' 'schema\tpackstar-preflight-sequence-v2'
    printf 'shardCount\t%s\n' "$SHARD_COUNT"
    printf 'totalBundles\t%s\n' "$TOTAL_BUNDLES"
    printf 'totalStates\t%s\n' "$TOTAL_STATES"
    printf 'replicatedBundleOrdinal\t%s\n' "$REPLICATED_BUNDLE_ORDINAL"
    printf 'scheduledBundles\t%s\n' "$TOTAL_SCHEDULED_BUNDLES"
    printf 'scheduledStates\t%s\n' "$TOTAL_ROWS"
    printf 'serialPredictedCaseHours\t%s\n' "$SERIAL_HOURS"
    printf 'multiNodeMakespanHours\t%s\n' "$MAKESPAN_HOURS"
    printf 'admissionMetric\t%s\n' "$ADMISSION_METRIC"
    printf 'admissionHours\t%s\n' "$ADMISSION_HOURS"
    printf 'slaHours\t%s\n' "$SLA_HOURS"
    printf 'admitted\t%s\n' "$ADMITTED"
    printf 'policyRows\t%s\n' "$POLICY_ROWS"
} > "$SUMMARY_TMP"
mv "$SUMMARY_TMP" "$REPORT_OUT"
mv "$DETAIL_TMP" "$DETAIL_OUT"

echo "packstar preflight aggregate: bundles=$TOTAL_BUNDLES states=$TOTAL_STATES serialHours=$SERIAL_HOURS makespanHours=$MAKESPAN_HOURS metric=$ADMISSION_METRIC admissionHours=$ADMISSION_HOURS slaHours=$SLA_HOURS admitted=$ADMITTED policies=$POLICY_ROWS"
if [ "$ADMITTED" != true ]; then
    exit 10
fi
