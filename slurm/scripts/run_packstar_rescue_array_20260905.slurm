#!/bin/bash
# Retry only the pfuncs selected by select_packstar_rescue_cases.py.
# Four A5000 GPUs and 48 CPU cores per task make 12 concurrent tasks fill the
# 48-GPU GrisMan A5000 pool while keeping the DP worker configuration identical
# within every task.  CCD remains on CPU here; the GPUs are used by branch-DP.
#SBATCH --job-name=ps_rescue_newproposal
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --gres=gpu:a5000:4
#SBATCH --time=14-00:00:00
#SBATCH --array=0-53%12
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=${BUILD_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/current_build}
CP_FILE="$BUILD_ROOT/test_classpath.txt"
SELECTION_ROOT=${SELECTION_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/selection_12505406}
CANDIDATES=${CANDIDATES:-$SELECTION_ROOT/candidates.tsv}
SEQUENCE_MAP_ROOT=${SEQUENCE_MAP_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/sequence_map_12505428}
EMAT_ROOT=${EMAT_ROOT:-/usr/xtmp/lz280/bench_comparison/results/emat_cache}
PDB_ROOT=${PDB_ROOT:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
RESULT_ROOT=${RESULT_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/newproposal_array_${SLURM_ARRAY_JOB_ID}}
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
THREADS=${SLURM_CPUS_PER_TASK:-48}
JAVA_HEAP_GB=${JAVA_HEAP_GB:-280}
DP_GPU_MAX_GPUS=${DP_GPU_MAX_GPUS:-4}

[ -s "$CP_FILE" ] || { echo "missing current-build classpath: $CP_FILE" >&2; exit 2; }
[ -s "$CANDIDATES" ] || { echo "missing candidate manifest: $CANDIDATES" >&2; exit 2; }
[ -d "$SEQUENCE_MAP_ROOT" ] || { echo "missing sequence map root: $SEQUENCE_MAP_ROOT" >&2; exit 2; }

candidate_line=$(awk -F '\t' -v n=$((SLURM_ARRAY_TASK_ID + 2)) \
    'NR == n { print; exit }' "$CANDIDATES")
[ -n "$candidate_line" ] || {
    echo "no candidate row for array task $SLURM_ARRAY_TASK_ID" >&2
    exit 2
}

case_id=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $1}')
design=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $2}')
pdb=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $3}')
mutable=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $4}')
flexible=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $5}')
state=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $6}')
sequence=$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $7}')

[ -n "$case_id" ] && [ -n "$design" ] && [ -n "$pdb" ] \
    && [ -n "$state" ] && [ -n "$sequence" ] || {
    echo "malformed candidate row: $candidate_line" >&2
    exit 2
}

sequence_map="$SEQUENCE_MAP_ROOT/${design}.tsv"
[ -s "$sequence_map" ] || { echo "missing sequence map: $sequence_map" >&2; exit 2; }
seq_index=$(awk -F '\t' -v s="$sequence" '$3 == s { print $2; exit }' "$sequence_map")
[ -n "$seq_index" ] || {
    echo "sequence not found in map: design=$design sequence=$sequence map=$sequence_map" >&2
    exit 2
}

OUT_DIR="$RESULT_ROOT/case_${case_id}_${design}_${state}_seq${seq_index}"
if [ -e "$OUT_DIR" ]; then
    echo "refusing to overwrite existing output: $OUT_DIR" >&2
    exit 2
fi
mkdir -p "$OUT_DIR/dp_mmap" "$OUT_DIR/emat_cache/$design" \
    "$OUT_DIR/adaptive_frequency_severity" /usr/xtmp/lz280/slurm_logs

SRC="$EMAT_ROOT/$design"
PDB_PATH="$PDB_ROOT/$pdb/$pdb.min.reduce.renum.pdb"
[ -d "$SRC" ] || { echo "missing emat source: $SRC" >&2; exit 2; }
[ -r "$PDB_PATH" ] || { echo "missing PDB: $PDB_PATH" >&2; exit 2; }
for state_name in complex ligand protein; do
    rigid="$SRC/export.$state_name.rigid.dat"
    minimizing="$SRC/export.$state_name.min.dat"
    [ -r "$rigid" ] && [ -r "$minimizing" ] || {
        echo "missing emat input: design=$design state=$state_name source=$SRC" >&2
        exit 2
    }
    ln -s "$rigid" "$OUT_DIR/emat_cache/$design/packstar.$state_name.rigid.dat"
    ln -s "$minimizing" "$OUT_DIR/emat_cache/$design/packstar.$state_name.minimizing.dat"
done

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
fi
gpu_inventory=unavailable
if command -v nvidia-smi >/dev/null 2>&1; then
    gpu_inventory=$(nvidia-smi --query-gpu=index,name,memory.total,uuid \
        --format=csv,noheader 2>/dev/null | paste -sd ';' - || true)
    [ -n "$gpu_inventory" ] || gpu_inventory=unavailable
fi

CP=$(<"$CP_FILE")
RESULT_TSV="$OUT_DIR/result.tsv"
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    "-Xmx${JAVA_HEAP_GB}g"
    -Xms4g
    -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$OUT_DIR/dp_mmap"
    -Dbranchdp.decomp.strategy=weighted_hicks
    -Dbranchdp.decomp.weightedHicks.restarts=24
    -Dbranchdp.decomp.weightedHicks.randomMoves=32
    -Dbranchdp.dp.parallel.threads="$THREADS"
    -Dpackstar.dp.parallel.threads="$THREADS"
    -Dbranchdp.dp.gpu=true
    -Dbranchdp.dp.gpu.failIfNoGpuPath=true
    -Dbranchdp.dp.gpu.failIfExceedsVram=true
    -Dbranchdp.dp.gpu.multiGpu=true
    -Dbranchdp.dp.gpu.maxGpus="$DP_GPU_MAX_GPUS"
    -Dbranchdp.dp.gpu.minMStatesPerGpu=4096
    -Dbranchdp.dp.gpu.childSliceMaxBytes=2147483648
    -Dbranchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dbranchdp.rootSplit=gpubytes
    -Dbranchdp.rootSplit.hostBudgetBytes=280000000000
    -Dbranchdp.rootSplit.gpuBudgetBytes=11000000000
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=0.5
    -Dpackstar.pac.frequencySeverity.jointMomentLearning=true
    -Dpackstar.pac.frequencySeverity.relativeBoundKcal=1
    -Dpackstar.pac.frequencySeverity.severityCap=20
    -Dpackstar.pac.frequencySeverity.severityPremiseId=conditional-relative-gauge-S0-20-not-externally-recalibrated
    -Dpackstar.pac.frequencySeverity.shrinkGrid=0:0,2:5,5:10,10:20
    -Dpackstar.pac.frequencySeverity.alphaGrid=1
    -Dpackstar.pac.frequencySeverity.folds=2
    -Dpackstar.pac.frequencySeverity.minShiftEssFraction=0.000000000001
    -Dpackstar.pac.frequencySeverity.discoveryMinShiftEssFraction=0.000000000001
    -Dpackstar.pac.frequencySeverity.maxRefits=8
    -Dpackstar.pac.frequencySeverity.discoverySamples=100
    -Dpackstar.pac.frequencySeverity.discoveryMaxSamples=400
    -Dpackstar.pac.frequencySeverity.validationSamples=400
    -Dpackstar.pac.frequencySeverity.tripleEta=true
    -Dpackstar.pac.frequencySeverity.tripleEtaScale=1
    -Dpackstar.pac.frequencySeverity.tripleEtaScaleGrid=1
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxAssignments=1000000
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxPositionTriples=3
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxFillEdges=3
    -Dpackstar.pac.frequencySeverity.tripleEtaMinCellContexts=1
    -Dpackstar.pac.frequencySeverity.tripleEtaPriorStrength=4
    -Dpackstar.pac.frequencySeverity.tripleEtaLocalCapKcal=2
    -Dpackstar.pac.frequencySeverity.tripleEtaResidualCapKcal=3
    -Dpackstar.pac.frequencySeverity.minTrainCount=5
    -Dpackstar.pac.frequencySeverity.maxUndertrainedAmplification=1.25
    -Dpackstar.pac.frequencySeverity.sizeSafety=0.9
    -Dpackstar.pac.frequencySeverity.severityTestAlpha=0.05
    -Dpackstar.pac.frequencySeverity.outputDir="$OUT_DIR/adaptive_frequency_severity"
    -Dpackstar.pac.samples=1000
    -Dpackstar.pac.trainSamples=500
    -Dpackstar.pac.pilotSamples=100
    -Dpackstar.pac.maxEstSamples=4000
    -Dpackstar.pac.monitorSamples=100
    -Dpackstar.pac.unreachableCap=400
    -Dpackstar.pac.nstarInflate=1.3
    -Dpackstar.pac.confidence=0.05
    -Dpackstar.pac.targetEpsilon=0.683
    -Dpackstar.pac.randomSeed=42
    -Dpackstar.pac.sampling.gpu=false
    -Dpackstar.pac.sampling.parallel=true
    -Dpackstar.pac.sampling.threads="$THREADS"
    -Dosprey.wmb.numGpus=0
    -Dosprey.wmb.streamsPerGpu=64
    -Dosprey.packstar.reduceMinimizations=true
    -Dosprey.packstar.correctionTightening=true
    -Dosprey.bench.method=packstar_pfunc
    -Dosprey.bench.epsilon=0.683
    -Dosprey.bench.numCPUs="$THREADS"
    -Dosprey.bench.outputDir="$OUT_DIR"
    -Dosprey.bench.designId="$design"
    -Dosprey.bench.pdbPath="$PDB_PATH"
    -Dosprey.bench.mutable="$mutable"
    -Dosprey.bench.flexible="$flexible"
    -Dosprey.packstarPfunc.state="$state"
    -Dosprey.packstarPfunc.seqIndex="$seq_index"
    -Dosprey.packstarPfunc.maxMut=1
    -Dosprey.packstarPfunc.outputTsv="$RESULT_TSV"
)

{
    printf 'schema\tpackstar-rescue-new-proposal-v1\n'
    printf 'arrayJob\t%s\n' "$SLURM_ARRAY_JOB_ID"
    printf 'arrayTask\t%s\n' "$SLURM_ARRAY_TASK_ID"
    printf 'caseId\t%s\n' "$case_id"
    printf 'design\t%s\n' "$design"
    printf 'pdb\t%s\n' "$pdb"
    printf 'state\t%s\n' "$state"
    printf 'sequenceIndex\t%s\n' "$seq_index"
    printf 'sequence\t%s\n' "$sequence"
    printf 'oldPackStarStatus\t%s\n' "$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $8}')"
    printf 'oldMarkStarStatus\t%s\n' "$(printf '%s\n' "$candidate_line" | awk -F '\t' '{print $9}')"
    printf 'proposalJointMomentLearning\ttrue\n'
    printf 'tripleMaximumPositionTriples\t3\n'
    printf 'tripleMaximumFillEdges\t3\n'
    printf 'tripleMaximumAssignments\t1000000\n'
    printf 'branchDpGpu\ttrue\n'
    printf 'branchDpGpuCount\t%s\n' "$DP_GPU_MAX_GPUS"
    printf 'ccdMode\tcpu\n'
    printf 'gres\tgpu:a5000:4\n'
    printf 'gpuInventory\t%s\n' "$gpu_inventory"
    printf 'startedAt\t%s\n' "$(date --iso-8601=seconds)"
} > "$OUT_DIR/run_manifest.tsv"

echo "rescue start job=$SLURM_JOB_ID array=$SLURM_ARRAY_JOB_ID task=$SLURM_ARRAY_TASK_ID case=$case_id design=$design state=$state seqIndex=$seq_index node=$SLURMD_NODENAME gpus=$gpu_inventory start=$(date --iso-8601=seconds)"

set +e
/usr/bin/time -o "$OUT_DIR/wall.time" \
    -f 'PACKSTAR_RESCUE_WALL elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" -cp "$CP" \
    edu.duke.cs.osprey.markstar.bench.GenericPDBBench \
    > >(tee "$OUT_DIR/run.log") 2>&1
java_status=${PIPESTATUS[0]}
set -e

printf 'javaStatus\t%s\n' "$java_status" >> "$OUT_DIR/run_manifest.tsv"
printf 'endedAt\t%s\n' "$(date --iso-8601=seconds)" >> "$OUT_DIR/run_manifest.tsv"
if [ -s "$RESULT_TSV" ]; then
    result_status=$(awk -F '\t' 'NR == 2 {print $5; exit}' "$RESULT_TSV")
    printf 'newPackStarStatus\t%s\n' "$result_status" >> "$OUT_DIR/run_manifest.tsv"
else
    printf 'newPackStarStatus\tNO_RESULT\n' >> "$OUT_DIR/run_manifest.tsv"
fi

echo "rescue end job=$SLURM_JOB_ID array=$SLURM_ARRAY_JOB_ID task=$SLURM_ARRAY_TASK_ID case=$case_id javaStatus=$java_status output=$OUT_DIR end=$(date --iso-8601=seconds)"
exit "$java_status"
