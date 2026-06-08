#!/bin/bash
# Submit compiled audit probe jobs from the 38-PDB benchmark spec.
set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=${OUTDIR:-/usr/xtmp/lz280/bench_comparison}
SPECS=${SPECS:-$OUTDIR/design_specs_prepped.csv}
PDBDIR=${PDBDIR:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
AUDIT_DIR=${AUDIT_DIR:-$OUTDIR/audit_leaves/gnn_s11_leafonly}
LOGDIR=${LOGDIR:-$OUTDIR/logs/compiled_audit_matrix}

# MODE=specs submits one job per design from design_specs_prepped.csv.
# MODE=leaves submits one job per design/state with existing S11 leaf CSVs.
MODE=${MODE:-specs}
STATE=${STATE:-Complex}
DESIGNS=${DESIGNS:-all}
LIMIT=${LIMIT:-0}
DRY_RUN=${DRY_RUN:-false}

COMPILE=${COMPILE:-true}
PREVIEW_CONFS=${PREVIEW_CONFS:-1}
COMPILE_THREADS=${COMPILE_THREADS:-1}
JAVA_XMX=${JAVA_XMX:-40g}
RUN_BUILD=${RUN_BUILD:-false}
CUDA_SWEEP=${CUDA_SWEEP:-false}
MAX_CONFS=${MAX_CONFS:-512}
WARMUP_CONFS=${WARMUP_CONFS:-64}
REPEATS=${REPEATS:-2}
STREAMS_LIST=${STREAMS_LIST:-"8,16,32,64,128"}
BATCH_SIZES=${BATCH_SIZES:-"256,512,1024,2048,4096"}
PRECISION=${PRECISION:-Float32}

PARTITION=${PARTITION:-grisman}
ACCOUNT=${ACCOUNT:-grisman}
if [ -z "${GRES+x}" ]; then
    if [ "$CUDA_SWEEP" = "true" ]; then
        GRES=gpu:a5000:1
    else
        GRES=
    fi
fi
NODELIST=${NODELIST:-fennario-01,fennario-02,fennario-03}
CPUS=${CPUS:-8}
MEM=${MEM:-64G}
TIME=${TIME:-02:00:00}
JOB_NAME_PREFIX=${JOB_NAME_PREFIX:-compiled_probe}

if [ "$DRY_RUN" != "true" ]; then
    mkdir -p "$LOGDIR"
fi

contains_design() {
    local did=$1
    if [ "$DESIGNS" = "all" ]; then
        return 0
    fi
    for wanted in $DESIGNS; do
        if [ "$wanted" = "$did" ]; then
            return 0
        fi
    done
    return 1
}

first_leaf_csv() {
    local did=$1
    local state=$2
    local dir="$AUDIT_DIR/$did/$state"
    if [ ! -d "$dir" ]; then
        return 0
    fi
    find "$dir" -maxdepth 1 -type f -name 'seq_*.csv' ! -name '*.summary.csv' 2>/dev/null | sort | head -1
}

submit_one() {
    local did=$1
    local pdb=$2
    local state=$3
    local mutable=$4
    local flexible=$5
    local input=$6

    local pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"
    if [ ! -f "$pdbpath" ]; then
        echo "SKIP $did/$state: missing PDB $pdbpath" >&2
        return
    fi
    if [ -z "$input" ]; then
        input="$AUDIT_DIR/$did/$state/no_leaf_csv_available.csv"
    fi

    local flags=(--parsable
        --partition="$PARTITION"
        --nodes=1
        --cpus-per-task="$CPUS"
        --mem="$MEM"
        --time="$TIME"
        --job-name="${JOB_NAME_PREFIX}_${did}_${state}")
    if [ -n "$ACCOUNT" ]; then
        flags+=(--account="$ACCOUNT")
    fi
    if [ -n "$GRES" ]; then
        flags+=(--gres="$GRES")
    fi
    if [ -n "$NODELIST" ]; then
        flags+=(--nodelist="$NODELIST")
    fi

    local wrap
    printf -v wrap 'cd /home/users/lz280/IdeaProjects/OSPREY3 && echo "compiled audit matrix did=%q state=%q host=$(hostname) date=$(date)" && DID=%q PDB=%q STATE=%q INPUT=%q MUT=%q FLEX=%q COMPILE=%q PREVIEW_CONFS=%q COMPILE_THREADS=%q JAVA_XMX=%q RUN_BUILD=%q CUDA_SWEEP=%q MAX_CONFS=%q WARMUP_CONFS=%q REPEATS=%q STREAMS_LIST=%q BATCH_SIZES=%q PRECISION=%q bash bench_comparison/probe_compiled_audit_mapping.sh' \
        "$did" "$state" "$did" "$pdbpath" "$state" "$input" "$mutable" "$flexible" "$COMPILE" "$PREVIEW_CONFS" "$COMPILE_THREADS" "$JAVA_XMX" "$RUN_BUILD" "$CUDA_SWEEP" "$MAX_CONFS" "$WARMUP_CONFS" "$REPEATS" "$STREAMS_LIST" "$BATCH_SIZES" "$PRECISION"

    if [ "$DRY_RUN" = "true" ]; then
        echo "DRY $did/$state input=$input"
    else
        local jid
        jid=$(sbatch "${flags[@]}" \
            --output="$LOGDIR/${did}_${state}_%j.out" \
            --error="$LOGDIR/${did}_${state}_%j.err" \
            --wrap "$wrap")
        echo "Submitted $did/$state: $jid -> $LOGDIR/${did}_${state}_${jid}.out"
    fi
}

submitted=0
while IFS=',' read -r did pdb desc nseq nres mutable flexible nconf trank; do
    if [ -z "$did" ] || [[ "$did" == \#* ]]; then
        continue
    fi
    if ! contains_design "$did"; then
        continue
    fi

    if [ "$MODE" = "leaves" ]; then
        if [ ! -d "$AUDIT_DIR/$did" ]; then
            continue
        fi
        while IFS= read -r state_dir; do
            state=$(basename "$state_dir")
            input=$(first_leaf_csv "$did" "$state")
            submit_one "$did" "$pdb" "$state" "$mutable" "$flexible" "$input"
            submitted=$((submitted + 1))
            if [ "$LIMIT" -gt 0 ] && [ "$submitted" -ge "$LIMIT" ]; then
                exit 0
            fi
        done < <(find "$AUDIT_DIR/$did" -mindepth 1 -maxdepth 1 -type d | sort)
    else
        input=$(first_leaf_csv "$did" "$STATE")
        submit_one "$did" "$pdb" "$STATE" "$mutable" "$flexible" "$input"
        submitted=$((submitted + 1))
        if [ "$LIMIT" -gt 0 ] && [ "$submitted" -ge "$LIMIT" ]; then
            exit 0
        fi
    fi
done < "$SPECS"

echo "Submitted $submitted compiled audit matrix job(s)"
