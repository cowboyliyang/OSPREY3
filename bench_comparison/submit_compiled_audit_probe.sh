#!/bin/bash
# Submit the compiled ConfSpace audit-mapping probe to Slurm.
set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=${OUTDIR:-/usr/xtmp/lz280/bench_comparison}
LOGDIR=${LOGDIR:-$OUTDIR/logs}

DID=${DID:-3u7y}
STATE=${STATE:-Complex}
PDB=${PDB:-/usr/xtmp/lz280/dance_bench/pdbs_prepped/$DID/${DID}.min.reduce.renum.pdb}
INPUT=${INPUT:-$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv}
MUT=${MUT:-"G384;G382"}
FLEX=${FLEX:-"G385;L587;L649"}
COMPILE=${COMPILE:-true}
PREVIEW_CONFS=${PREVIEW_CONFS:-1}
COMPILE_THREADS=${COMPILE_THREADS:-1}
JAVA_XMX=${JAVA_XMX:-40g}
RUN_BUILD=${RUN_BUILD:-true}
CUDA_SWEEP=${CUDA_SWEEP:-false}
MAX_CONFS=${MAX_CONFS:-512}
WARMUP_CONFS=${WARMUP_CONFS:-64}
REPEATS=${REPEATS:-2}
STREAMS_LIST=${STREAMS_LIST:-"8,16,32,64,128"}
BATCH_SIZES=${BATCH_SIZES:-"256,512,1024,2048,4096"}
PRECISION=${PRECISION:-Float32}

PARTITION=${PARTITION:-grisman}
ACCOUNT=${ACCOUNT:-grisman}
GRES=${GRES:-gpu:a5000:1}
NODELIST=${NODELIST:-fennario-01,fennario-02,fennario-03}
CPUS=${CPUS:-8}
MEM=${MEM:-64G}
TIME=${TIME:-02:00:00}
JOB_NAME=${JOB_NAME:-compiled_probe}

mkdir -p "$LOGDIR"

SBATCH_FLAGS=(--parsable
    --partition="$PARTITION"
    --nodes=1
    --cpus-per-task="$CPUS"
    --mem="$MEM"
    --time="$TIME"
    --job-name="$JOB_NAME")
if [ -n "$ACCOUNT" ]; then
    SBATCH_FLAGS+=(--account="$ACCOUNT")
fi
if [ -n "$GRES" ]; then
    SBATCH_FLAGS+=(--gres="$GRES")
fi
if [ -n "$NODELIST" ]; then
    SBATCH_FLAGS+=(--nodelist="$NODELIST")
fi

printf -v WRAP 'cd /home/users/lz280/IdeaProjects/OSPREY3 && echo "compiled audit probe on $(hostname) at $(date)" && DID=%q PDB=%q STATE=%q INPUT=%q MUT=%q FLEX=%q COMPILE=%q PREVIEW_CONFS=%q COMPILE_THREADS=%q JAVA_XMX=%q RUN_BUILD=%q CUDA_SWEEP=%q MAX_CONFS=%q WARMUP_CONFS=%q REPEATS=%q STREAMS_LIST=%q BATCH_SIZES=%q PRECISION=%q bash bench_comparison/probe_compiled_audit_mapping.sh' \
    "$DID" "$PDB" "$STATE" "$INPUT" "$MUT" "$FLEX" "$COMPILE" "$PREVIEW_CONFS" "$COMPILE_THREADS" "$JAVA_XMX" "$RUN_BUILD" "$CUDA_SWEEP" "$MAX_CONFS" "$WARMUP_CONFS" "$REPEATS" "$STREAMS_LIST" "$BATCH_SIZES" "$PRECISION"

JID=$(sbatch "${SBATCH_FLAGS[@]}" \
    --output="$LOGDIR/compiled_audit_probe_%j.out" \
    --error="$LOGDIR/compiled_audit_probe_%j.err" \
    --wrap "$WRAP")

echo "Submitted compiled audit probe: $JID -> $LOGDIR/compiled_audit_probe_${JID}.out"
