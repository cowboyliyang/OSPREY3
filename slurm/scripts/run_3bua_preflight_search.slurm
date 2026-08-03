#!/bin/bash
# 3bua-only PACK* admission preflight with the improved multi-round + parallel
# branchwidth/root search.  Allocation-free, CPU-only (A5000 profile plugged in).
#SBATCH --job-name=3bua_pf_search
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

# --- A5000 profile (numbers only; no GPU is used by the dry-run preflight) ---
HW_LABEL=A5000_8x
GPU_RATE_PER_GPU=614000000
GPU_COUNT=8
OOC_RATE=170000000
GPU_BUDGET_BYTES=23085449216

# --- admission / search budget ---
CASE_SLA_HOURS=336
SOFT_STATE_HOURS=24
SAFETY_FACTOR=1.5
INITIAL_EXACT_MILLIS=300000     # 5 min initial exact-improve per state
FINAL_EXACT_MILLIS=600000       # 10 min per deepening attempt
FINAL_MAX_STATES=79             # allow every unique state to be deepened
MAX_ROUNDS=8                    # multi-round escalating branchwidth drop
PREVIEW_THREADS=96              # parallel per-state preview

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
CACHE_ROOT=/usr/xtmp/lz280/bench_comparison/results/emat_cache
RUN_ROOT=/tmp/packstar-3bua-pf-${SLURM_JOB_ID}
HOST_BUDGET_BYTES=450971566080
POLICY_DIR=/usr/xtmp/lz280/packstar_production/preflight_policies
BUDGET=${BUDGET:-0.5}
POLICY_OUT=$POLICY_DIR/3bua_rb${BUDGET}_A5000.tsv

cd "$REPO"
mkdir -p "$RUN_ROOT/emat_cache" "$POLICY_DIR"
cp -a "$CACHE_ROOT/3bua_fennario" "$RUN_ROOT/emat_cache/"

CP=$(<"$CP_FILE")
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx320g
    -Xms8g
    -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$RUN_ROOT/dp_mmap"
    -Dpackstar.decomp.strategy=adaptive
    -Dpackstar.decomp.weightedHicks.restarts=48
    -Dpackstar.decomp.weightedHicks.randomMoves=32
    -Dpackstar.decomp.exactImprove.maxDrop=3
    -Dpackstar.decomp.exactImprove.maxMillis="$INITIAL_EXACT_MILLIS"
    -Dpackstar.rootSplit=predicted
    -Dpackstar.rootSplit.hostBudgetBytes="$HOST_BUDGET_BYTES"
    -Dpackstar.rootSplit.gpuBudgetBytes="$GPU_BUDGET_BYTES"
    -Dpackstar.dp.gpu=true
    -Dpackstar.dp.gpu.childSliceMaxBytes=2147483648
    -Dpackstar.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dpackstar.admission.gpuWorkPerSecondPerGpu="$GPU_RATE_PER_GPU"
    -Dpackstar.admission.gpuCount="$GPU_COUNT"
    -Dpackstar.admission.oocBytesPerSecond="$OOC_RATE"
    -Dpackstar.admission.safetyFactor="$SAFETY_FACTOR"
    -Dpackstar.admission.softStateHours="$SOFT_STATE_HOURS"
    -Dpackstar.admission.caseSlaHours="$CASE_SLA_HOURS"
    -Dpackstar.admission.finalMaxStates="$FINAL_MAX_STATES"
    -Dpackstar.admission.finalExactMaxMillis="$FINAL_EXACT_MILLIS"
    -Dpackstar.admission.maxRounds="$MAX_ROUNDS"
    -Dpackstar.admission.previewThreads="$PREVIEW_THREADS"
    -Dpackstar.admission.policyOut="$POLICY_OUT"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget="$BUDGET"
    -Dosprey.bench.method=packstar
    -Dosprey.bench.packstarPreflightOnly=true
    -Dosprey.bench.numCPUs=96
    -Dosprey.bench.outputDir="$RUN_ROOT"
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

echo "preflight_job=$SLURM_JOB_ID hardware=$HW_LABEL gpuRatePerGpu=$GPU_RATE_PER_GPU gpuCount=$GPU_COUNT oocRate=$OOC_RATE gpuBudgetBytes=$GPU_BUDGET_BYTES safetyFactor=$SAFETY_FACTOR softStateHours=$SOFT_STATE_HOURS caseSlaHours=$CASE_SLA_HOURS residualBudget=$BUDGET maxRounds=$MAX_ROUNDS previewThreads=$PREVIEW_THREADS finalMaxStates=$FINAL_MAX_STATES finalExactMillis=$FINAL_EXACT_MILLIS start=$(date --iso-8601=seconds)"
echo "git_head=$(git rev-parse HEAD)"
sha256sum \
    build/classes/java/main/edu/duke/cs/osprey/packstar/PackStarAdmissionDecision.class \
    build/classes/java/main/edu/duke/cs/osprey/packstar/PackStarCasePreflight.class

status=0
echo "BEGIN case=3bua hardware=$HW_LABEL time=$(date --iso-8601=seconds)"
if /usr/bin/time -f "PREFLIGHT_WALL hardware=$HW_LABEL design=3bua elapsed=%e maxRssKiB=%M" \
    "$JAVA" "${JARGS[@]}" \
    -Dosprey.bench.designId=3bua_fennario \
    -Dosprey.bench.pdbPath="$PDB_ROOT/3bua/3bua.min.reduce.renum.pdb" \
    '-Dosprey.bench.mutable=C447;C450' \
    '-Dosprey.bench.flexible=D638;D641;D642;D645;D623;D619;D616;D620;C417;C414;C410;C421;C452;C443;C444;C446' \
    -cp "$CP" "$MAIN"; then
    echo "END case=3bua hardware=$HW_LABEL status=admitted time=$(date --iso-8601=seconds)"
else
    rc=$?
    status=1
    echo "END case=3bua hardware=$HW_LABEL status=rejected-or-failed rc=$rc time=$(date --iso-8601=seconds)"
fi
echo "preflight_job=$SLURM_JOB_ID hardware=$HW_LABEL end=$(date --iso-8601=seconds) status=$status"
exit "$status"
