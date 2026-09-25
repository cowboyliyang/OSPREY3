#!/bin/bash
#SBATCH --job-name=packstar-newprop-38
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --array=0-37%38
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --gres=gpu:4
#SBATCH --time=14-00:00:00
#SBATCH --output=/usr/xtmp/lz280/packstar_adaptive_frequency_severity/slurm_logs/newproposal_38_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/packstar_adaptive_frequency_severity/slurm_logs/newproposal_38_%A_%a.err

set -euo pipefail

REPO="/home/users/lz280/OSPREY3-fresh-packstar"
BUILD_ROOT="/usr/xtmp/lz280/packstar_rescue_20260905/current_build"
CP_FILE="$BUILD_ROOT/test_classpath.txt"
JAVA="/home/users/lz280/java/jdk-17.0.2+8/bin/java"
CSV="/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv"
PDB_ROOT="/usr/xtmp/lz280/dance_bench/pdbs_prepped"
EMAT_ROOT="/usr/xtmp/lz280/bench_comparison/results/emat_cache"
PRIOR_ROOT="/usr/xtmp/lz280/packstar_adaptive_frequency_severity/decomposition_cost_grisman_multigpu4_38_A12494208"

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:?missing SLURM_ARRAY_JOB_ID}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"
RESULT_ROOT="${RESULT_ROOT:-/usr/xtmp/lz280/packstar_adaptive_frequency_severity/newproposal_joint_moment_dedup_38_A${ARRAY_JOB_ID}}"

mkdir -p "$RESULT_ROOT/slurm_logs"

# The design CSV starts with two comment/header lines.  Map array indices to
# the 38 data rows in the same order as the previous full-system run.
ROW="$(awk -F, -v n="$((TASK_ID + 1))" '
    NR > 2 && $1 !~ /^#/ && NF > 0 {
        count++
        if (count == n) {
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

OUT_DIR="$RESULT_ROOT/runs/${DESIGN_ID}_A${ARRAY_JOB_ID}_T${TASK_ID}"
RANK_DIR="$RESULT_ROOT/rankings"
mkdir -p "$OUT_DIR/emat_cache/$DESIGN_ID" "$OUT_DIR/dp_mmap" "$RANK_DIR"

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
    ln -sfn "$EMAT_ROOT/$DESIGN_ID/export.$state.rigid.dat" \
        "$OUT_DIR/emat_cache/$DESIGN_ID/packstar.$state.rigid.dat"
    ln -sfn "$EMAT_ROOT/$DESIGN_ID/export.$state.min.dat" \
        "$OUT_DIR/emat_cache/$DESIGN_ID/packstar.$state.minimizing.dat"
done

source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true

GPU_NAME="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1 || true)"
GPU_COUNT="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ' || true)"
if [[ -z "$GPU_NAME" || "$GPU_COUNT" -lt 1 ]]; then
    echo "No CUDA GPU visible on allocated node" >&2
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

GIT_HEAD="$(git -C "$REPO" rev-parse HEAD 2>/dev/null || true)"
BUILD_MANIFEST="$BUILD_ROOT/build_manifest.tsv"
{
    printf 'schema\tpackstar-newproposal-38-v1\n'
    printf 'designId\t%s\n' "$DESIGN_ID"
    printf 'pdbId\t%s\n' "$PDB_ID"
    printf 'arrayJobId\t%s\n' "$ARRAY_JOB_ID"
    printf 'arrayTaskId\t%s\n' "$TASK_ID"
    printf 'gitHead\t%s\n' "$GIT_HEAD"
    printf 'buildRoot\t%s\n' "$BUILD_ROOT"
    printf 'buildManifest\t%s\n' "$BUILD_MANIFEST"
    printf 'prior38Root\t%s\n' "$PRIOR_ROOT"
    printf 'proposalLearning\tjoint-pair-fit-and-triple-second-moment\n'
    printf 'jointMomentLearning\ttrue\n'
    printf 'conformationDedup\texact-automatic\n'
    printf 'tripleMaximumPositionTriples\t3\n'
    printf 'tripleMaximumFillEdges\t3\n'
    printf 'tripleMaximumAssignments\t1000000\n'
    printf 'ccdMode\tcpu\n'
    printf 'branchDpGpu\ttrue\n'
    printf 'branchDpGpuCount\t%s\n' "$GPU_COUNT"
    printf 'gpuName\t%s\n' "$GPU_NAME"
    printf 'gpuBudgetBytes\t%s\n' "$GPU_BUDGET_BYTES"
    printf 'gres\tgpu:4\n'
    printf 'targetEpsilon\t0.683\n'
} > "$OUT_DIR/run_manifest.tsv"

JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx280g
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
    -Dpackstar.pac.randomSeed=42
    -Dpackstar.pac.sampling.gpu=false
    -Dpackstar.pac.sampling.parallel=true
    -Dpackstar.pac.sampling.threads=48
    -Dosprey.wmb.numGpus=0
    -Dosprey.wmb.streamsPerGpu=64
    -Dosprey.packstar.reduceMinimizations=true
    -Dosprey.packstar.correctionTightening=true
    -Dosprey.bench.method=packstar
    -Dosprey.bench.epsilon=0.683
    -Dosprey.bench.numCPUs=48
    "-Dosprey.bench.outputDir=$OUT_DIR"
    "-Dosprey.bench.designId=$DESIGN_ID"
    "-Dosprey.bench.pdbPath=$PDB_PATH"
    "-Dosprey.bench.mutable=$MUTABLE"
    "-Dosprey.bench.flexible=$FLEXIBLE"
)

echo "PACKSTAR-38 start design=$DESIGN_ID pdb=$PDB_ID task=$TASK_ID gpu=$GPU_NAME count=$GPU_COUNT"
echo "PACKSTAR-38 resultRoot=$RESULT_ROOT"
echo "PACKSTAR-38 proposal=joint-moment dedup=exact tripleMax=3"

set +e
/usr/bin/time -o "$OUT_DIR/wall.time" -f 'PACKSTAR_38_WALL elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" -cp "$(cat "$CP_FILE")" \
    edu.duke.cs.osprey.markstar.bench.GenericPDBBench \
    > >(tee "$OUT_DIR/run.log") 2>&1
JAVA_STATUS=${PIPESTATUS[0]}
set -e

RESULT_CSV="$OUT_DIR/${DESIGN_ID}_packstar.csv"
if [[ -s "$RESULT_CSV" ]]; then
    cp "$RESULT_CSV" "$RANK_DIR/${DESIGN_ID}_packstar.csv"
    RESULT_PRESENT=true
else
    RESULT_PRESENT=false
fi

DEDUP_LINES="$(grep -c 'exact conformation dedup' "$OUT_DIR/run.log" 2>/dev/null || true)"
JOINT_PAIR_LINES="$(grep -c 'joint-pair-fit' "$OUT_DIR/run.log" 2>/dev/null || true)"
TRIPLE_M2_LINES="$(grep -c 'triple-m2-fit' "$OUT_DIR/run.log" 2>/dev/null || true)"
if [[ -s "$RESULT_CSV" ]]; then
    ESTIMATED_ROWS="$(awk -F, 'NR > 1 && $0 ~ /,Estimated,/ { n++ } END { print n + 0 }' "$RESULT_CSV")"
    ABORTED_ROWS="$(awk -F, 'NR > 1 && $0 ~ /,Aborted,/ { n++ } END { print n + 0 }' "$RESULT_CSV")"
    OTHER_ROWS="$(awk -F, 'NR > 1 { if ($0 !~ /,Estimated,/ && $0 !~ /,Aborted,/) n++ } END { print n + 0 }' "$RESULT_CSV")"
else
    ESTIMATED_ROWS=0
    ABORTED_ROWS=0
    OTHER_ROWS=0
fi

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

# Mmap files are generated intermediates and can be very large.  Keep all
# reports/logs, but remove only this task's explicitly generated workspace.
if [[ "${KEEP_MMAP:-false}" != true ]]; then
    rm -rf -- "$OUT_DIR/dp_mmap"
fi

echo "PACKSTAR-38 done design=$DESIGN_ID javaStatus=$JAVA_STATUS estimated=$ESTIMATED_ROWS aborted=$ABORTED_ROWS"
