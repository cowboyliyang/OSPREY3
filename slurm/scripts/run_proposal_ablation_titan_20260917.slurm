#!/bin/bash
#SBATCH --job-name=ps_proposal_titan
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --exclude=jerry3
#SBATCH --array=0-23%7
#SBATCH --time=1-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
TASK_INDEX=${SLURM_ARRAY_TASK_ID:-0}
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}
SYSTEMS=(1gwc 2p4a 2xxm 3k3q 2xgy 3ma2 4wyu 5a6y)
SELECTED_SYSTEM=${CASE_OVERRIDE:-${SYSTEMS[$((TASK_INDEX / 3))]}}
SEED=${SEED_OVERRIDE:-$((42 + TASK_INDEX % 3))}
RESULT_ROOT=${RESULT_ROOT_OVERRIDE:-/usr/xtmp/lz280/packstar_proposal_ablation_20260917/A${ARRAY_JOB_ID}}
RUNNER_SOURCE=${BASH_SOURCE[0]}
test -f "${BUILD_ROOT:?}/READY"
# Each case runs every arm sequentially on the same allocated node/GPUs.
# Rotate arm order by case; seed, caches, CPU CCD, and all PAC controls are fixed.
run_arm() (
set -euo pipefail

REPO="${BUILD_ROOT:?test build root required}/source"
BUILD_ROOT="${BUILD_ROOT:?}"
cd "$REPO"
CP_FILE="$BUILD_ROOT/test_classpath.txt"
JAVA="/home/users/lz280/java/jdk-17.0.2+8/bin/java"
CSV="/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv"
PDB_ROOT="/usr/xtmp/lz280/dance_bench/pdbs_prepped"
EMAT_ROOT="/usr/xtmp/lz280/bench_comparison/results/emat_cache"
PRIOR_ROOT="/usr/xtmp/lz280/packstar_adaptive_frequency_severity/decomposition_cost_grisman_multigpu4_38_A12494208"

TASK_ID=$TASK_INDEX

mkdir -p "$RESULT_ROOT/slurm_logs"

# The design CSV starts with two comment/header lines.  Map array indices to
# the 38 data rows in the same order as the previous full-system run.
ROW="$(awk -F, -v selected="$SELECTED_SYSTEM" '
    NR > 2 && $1 !~ /^#/ && NF > 0 {
        if ($1 == selected) {
            print
            exit
        }
    }
' "$CSV")"

if [[ -z "$ROW" ]]; then
    echo "No design row for array task $TASK_ID" >&2
    exit 2
fi

IFS=',' read -r DESIGN_ID PDB_ID _ _ _ MUTABLE FLEXIBLE _ <<< "$ROW"
if [[ -z "$DESIGN_ID" || -z "$PDB_ID" || -z "$MUTABLE" || -z "$FLEXIBLE" ]]; then
    echo "Malformed design row: $ROW" >&2
    exit 2
fi

OUT_DIR="$RESULT_ROOT/$DESIGN_ID/seed_$SEED/$ARM"
RANK_DIR="$OUT_DIR/rankings"
mkdir -p "$OUT_DIR/emat_cache/$DESIGN_ID" "$OUT_DIR/dp_mmap" "$RANK_DIR" "$OUT_DIR/tmp"
export TMPDIR="$OUT_DIR/tmp"
export CUDA_CACHE_PATH="$OUT_DIR/cuda_cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "$RUNNER_SOURCE" "$OUT_DIR/runner.slurm"
printf '%s\n' "$ROW" > "$OUT_DIR/design_row.csv"

PDB_PATH="$PDB_ROOT/$PDB_ID/$PDB_ID.min.reduce.renum.pdb"
if [[ ! -s "$PDB_PATH" ]]; then
    echo "Missing PDB: $PDB_PATH" >&2
    exit 3
fi
if [[ ! -s "$CP_FILE" ]]; then
    echo "Missing current build classpath: $CP_FILE" >&2
    exit 3
fi

# PackStar expects these six names in the per-run cache directory.
for state in complex ligand protein; do
    cp --reflink=auto "$EMAT_ROOT/$DESIGN_ID/export.$state.rigid.dat" \
        "$OUT_DIR/emat_cache/$DESIGN_ID/packstar.$state.rigid.dat"
    cp --reflink=auto "$EMAT_ROOT/$DESIGN_ID/export.$state.min.dat" \
        "$OUT_DIR/emat_cache/$DESIGN_ID/packstar.$state.minimizing.dat"
done

source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true

GPU_NAME="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1 || true)"
GPU_COUNT="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ' || true)"
if [[ "$GPU_NAME" != *"TITAN V"* || "$GPU_COUNT" -ne 4 || "${SLURM_CPUS_PER_TASK:-0}" -ne 48 ]]; then
    echo "Expected 48 CPUs and four TITAN V GPUs" >&2
    exit 4
fi

# Match the previous run's hardware-specific root-split budget while keeping
# the job usable if GrisMan assigns another supported GPU model.
GPU_BUDGET_BYTES=11000000000
GPU_RATE=614000000
GPU_OOC_RATE=170000000
if [[ "$GPU_NAME" == *"A5000"* ]]; then
    GPU_BUDGET_BYTES=23085449216
    GPU_RATE=614000000
    GPU_OOC_RATE=170000000
elif [[ "$GPU_NAME" == *"TITAN V"* || "$GPU_NAME" == *"Titan V"* ]]; then
    GPU_BUDGET_BYTES=10737418240
    GPU_RATE=1697000000
    GPU_OOC_RATE=190000000
fi

GIT_HEAD="source-sha256-manifest"
BUILD_MANIFEST="$BUILD_ROOT/source.sha256"
nvidia-smi --query-gpu=index,name,uuid,memory.total,driver_version --format=csv > "$OUT_DIR/gpus.csv"
lscpu > "$OUT_DIR/lscpu.txt"
scontrol show job "$SLURM_JOB_ID" > "$OUT_DIR/slurm_job.txt"
sha256sum "$PDB_PATH" "$OUT_DIR/emat_cache/$DESIGN_ID/"*.dat > "$OUT_DIR/input.sha256"
{
    printf 'schema\tpackstar-proposal-ablation-v1\n'
    printf 'seed\t%s\n' "$SEED"
    printf 'node\t%s\n' "$SLURMD_NODENAME"
    printf 'precision\tFP64\n'
    printf 'proposalLearningEnabled\t%s\n' "$LEARNING"
    printf 'tripleEtaEnabled\t%s\n' "$TRIPLES"
    printf 'selector\t%s\n' "$STRATEGY"
    printf 'ematMode\tprivate-copy-historical-cache-excluded-from-wall\n'
    printf 'cpuThreads\t48\n'
    printf 'hostBudgetBytes\t280000000000\n'
    printf 'arm\t%s\n' "$ARM"
    printf 'armOrder\t%s\n' "$ARM_ORDER"
    printf 'scope\tcomplex-seq%s\n' "${SEQ_OVERRIDE:-0}"
    printf 'designId\t%s\n' "$DESIGN_ID"
    printf 'pdbId\t%s\n' "$PDB_ID"
    printf 'arrayJobId\t%s\n' "$ARRAY_JOB_ID"
    printf 'arrayTaskId\t%s\n' "$TASK_ID"
    printf 'gitHead\t%s\n' "$GIT_HEAD"
    printf 'buildRoot\t%s\n' "$BUILD_ROOT"
    printf 'buildManifest\t%s\n' "$BUILD_MANIFEST"
    printf 'prior38Root\t%s\n' "$PRIOR_ROOT"
    printf 'conformationDedup\texact-automatic\n'
    printf 'tripleMaximumPositionTriples\t3\n'
    printf 'tripleMaximumFillEdges\t3\n'
    printf 'tripleMaximumAssignments\t1000000\n'
    printf 'ccdMode\tcpu\n'
    printf 'branchDpGpu\ttrue\n'
    printf 'branchDpGpuCount\t%s\n' "$GPU_COUNT"
    printf 'gpuName\t%s\n' "$GPU_NAME"
    printf 'gpuBudgetBytes\t%s\n' "$GPU_BUDGET_BYTES"
    printf 'gres\tgpu:titan_v:4\n'
    printf 'targetEpsilon\t0.683\n'
} > "$OUT_DIR/run_manifest.tsv"

JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx280g
    -Xms4g
    -XX:-UseSuperWord
    "-Djava.io.tmpdir=$TMPDIR"

    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$OUT_DIR/dp_mmap"
    -Dbranchdp.decomp.strategy=weighted_hicks
    -Dbranchdp.decomp.weightedHicks.restarts=24
    -Dbranchdp.decomp.weightedHicks.randomMoves=32
    -Dbranchdp.dp.parallel.threads=48
    -Dpackstar.dp.parallel.threads=48
    -Dbranchdp.dp.gpu=true
    -Dbranchdp.dp.gpu.failIfNoGpuPath=true
    -Dbranchdp.dp.gpu.failIfExceedsVram=true
    -Dbranchdp.dp.gpu.multiGpu=true
    -Dbranchdp.dp.gpu.maxGpus=4
    -Dbranchdp.dp.gpu.minMStatesPerGpu=4096
    -Dbranchdp.dp.gpu.childSliceMaxBytes=2147483648
    -Dbranchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dbranchdp.rootSplit=gpubytes
    -Dbranchdp.rootSplit.hostBudgetBytes=280000000000
    "-Dbranchdp.rootSplit.gpuBudgetBytes=$GPU_BUDGET_BYTES"
    "-Dpackstar.admission.gpuWorkPerSecondPerGpu=$GPU_RATE"
    "-Dpackstar.admission.gpuOutOfCoreWorkPerSecondPerGpu=$GPU_OOC_RATE"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=1

    "-Dpackstar.pac.frequencySeverity.proposalLearning=$LEARNING"
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
    "-Dpackstar.pac.frequencySeverity.tripleEta=$TRIPLES"
    "-Dpackstar.pac.frequencySeverity.tripleEtaSelectionStrategy=$STRATEGY"
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxWork=1000000000000
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxMStates=1000000000
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxTableBytes=64GiB
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxTotalTableBytes=256GiB
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxHostBytes=280000000000
    -Dpackstar.pac.frequencySeverity.tripleEtaMaxFileBytes=512GiB
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
    "-Dpackstar.pac.frequencySeverity.outputDir=$OUT_DIR/adaptive_frequency_severity"

    -Dpackstar.pac.samples=1000
    -Dpackstar.pac.trainSamples=500
    -Dpackstar.pac.pilotSamples=100
    -Dpackstar.pac.maxEstSamples=4000
    -Dpackstar.pac.monitorSamples=100
    -Dpackstar.pac.unreachableCap=400
    -Dpackstar.pac.nstarInflate=1.3
    -Dpackstar.pac.confidence=0.05
    -Dpackstar.pac.targetEpsilon=0.683
    "-Dpackstar.pac.randomSeed=$SEED"
    -Dpackstar.pac.sampling.gpu=false
    -Dpackstar.pac.sampling.parallel=true
    -Dpackstar.pac.sampling.threads=48
    -Dosprey.wmb.numGpus=0
    -Dosprey.wmb.streamsPerGpu=64
    -Dosprey.packstar.reduceMinimizations=true
    -Dosprey.packstar.correctionTightening=true
    -Dosprey.bench.method=packstar_pfunc
    -Dosprey.packstarPfunc.state=complex
    "-Dosprey.packstarPfunc.seqIndex=${SEQ_OVERRIDE:-0}"
    -Dpackstar.pac.sampling.timing=true
    -Dpackstar.pac.dp.timing=true
    -Dosprey.bench.epsilon=0.683
    -Dosprey.bench.numCPUs=48
    "-Dosprey.bench.outputDir=$OUT_DIR"
    "-Dosprey.bench.designId=$DESIGN_ID"
    "-Dosprey.bench.pdbPath=$PDB_PATH"
    "-Dosprey.bench.mutable=$MUTABLE"
    "-Dosprey.bench.flexible=$FLEXIBLE"
)
if [[ -n "${FIT_THREADS:-}" ]]; then
    JARGS+=("-Dpackstar.pac.frequencySeverity.tripleFitThreads=$FIT_THREADS")
fi

printf '%q ' "$JAVA" "${JARGS[@]}" -cp "$(cat "$CP_FILE")" edu.duke.cs.osprey.markstar.bench.GenericPDBBench > "$OUT_DIR/command.sh"
printf '\n' >> "$OUT_DIR/command.sh"
echo "PROPOSAL_ABLATION start design=$DESIGN_ID seed=$SEED arm=$ARM gpu=$GPU_NAME count=$GPU_COUNT"
echo "PACKSTAR-38 resultRoot=$RESULT_ROOT"

set +e
/usr/bin/time -o "$OUT_DIR/wall.time" -f 'PACKSTAR_38_WALL elapsed=%e maxRssKiB=%M' \
    timeout --signal=TERM --kill-after=60s 5h "$JAVA" "${JARGS[@]}" -cp "$(cat "$CP_FILE")" \
    edu.duke.cs.osprey.markstar.bench.GenericPDBBench \
    > >(tee "$OUT_DIR/run.log") 2>&1
JAVA_STATUS=$?
set -e

RESULT_CSV="$OUT_DIR/${DESIGN_ID}_packstar_pfunc.tsv"
if [[ -s "$RESULT_CSV" ]]; then
    cp "$RESULT_CSV" "$RANK_DIR/${DESIGN_ID}_packstar_pfunc.tsv"
    RESULT_PRESENT=true
else
    RESULT_PRESENT=false
fi

DEDUP_LINES="$(grep -c 'exact conformation dedup' "$OUT_DIR/run.log" 2>/dev/null || true)"
JOINT_PAIR_LINES="$(grep -c 'joint-pair-fit' "$OUT_DIR/run.log" 2>/dev/null || true)"
TRIPLE_M2_LINES="$(grep -c 'triple-m2-fit' "$OUT_DIR/run.log" 2>/dev/null || true)"
ESTIMATED_ROWS=NA
ABORTED_ROWS=NA
OTHER_ROWS=NA

{
    printf 'javaStatus\t%s\n' "$JAVA_STATUS"
    printf 'resultPresent\t%s\n' "$RESULT_PRESENT"
    printf 'estimatedRows\t%s\n' "$ESTIMATED_ROWS"
    printf 'abortedRows\t%s\n' "$ABORTED_ROWS"
    printf 'otherRows\t%s\n' "$OTHER_ROWS"
    printf 'dedupLogLines\t%s\n' "$DEDUP_LINES"
    printf 'jointPairFitLogLines\t%s\n' "$JOINT_PAIR_LINES"
    printf 'tripleM2FitLogLines\t%s\n' "$TRIPLE_M2_LINES"
    printf 'endedAt\t%s\n' "$(date -Is)"
} >> "$OUT_DIR/run_manifest.tsv"

# The output can legitimately contain Aborted sequence rows; a hard Java
# failure or a missing full-system CSV is the task failure condition.
if [[ "$JAVA_STATUS" -ne 0 ]]; then
    exit "$JAVA_STATUS"
fi
if [[ "$RESULT_PRESENT" != true ]]; then
    exit 5
fi

echo "PACKSTAR-38 done design=$DESIGN_ID javaStatus=$JAVA_STATUS estimated=$ESTIMATED_ROWS aborted=$ABORTED_ROWS"
rm -rf -- "$OUT_DIR/dp_mmap"

)
ARMS=(no-learning pair-only decomposition-cost budget-forward)
FAILED=0
for ARM_ORDER in 0 1 2 3; do
    ARM=${ARMS[$(((TASK_INDEX + ARM_ORDER) % 4))]}
    if [[ -n "${ARM_FILTER:-}" && "$ARM" != "$ARM_FILTER" ]]; then continue; fi
    LEARNING=true
    TRIPLES=true
    STRATEGY=decomposition-cost
    if [[ "$ARM" == no-learning ]]; then LEARNING=false; TRIPLES=false; fi
    if [[ "$ARM" == pair-only ]]; then TRIPLES=false; fi
    if [[ "$ARM" == budget-forward ]]; then STRATEGY=budget-forward; fi
    set +e
    ( set -e; run_arm )
    STATUS=$?
    set -e
    if [[ "$STATUS" != 0 ]]; then FAILED=1; fi
    echo "COMPARISON_ARM_EXIT case=$SELECTED_SYSTEM seed=$SEED arm=$ARM status=$STATUS"
done
exit "$FAILED"
