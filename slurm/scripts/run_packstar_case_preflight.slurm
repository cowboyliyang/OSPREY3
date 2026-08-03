#!/bin/bash
# Required exports: HW_LABEL, GPU_RATE_PER_GPU, GPU_COUNT, OOC_RATE,
# GPU_BUDGET_BYTES.  Submit after the matching calibration has completed.
#SBATCH --job-name=packstar_preflight
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00

set -euo pipefail

: "${HW_LABEL:?HW_LABEL is required}"
: "${GPU_RATE_PER_GPU:?GPU_RATE_PER_GPU is required}"
: "${GPU_COUNT:?GPU_COUNT is required}"
: "${OOC_RATE:?OOC_RATE is required}"
: "${GPU_BUDGET_BYTES:?GPU_BUDGET_BYTES is required}"

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
CACHE_ROOT=/usr/xtmp/lz280/bench_comparison/results/emat_cache
RUN_ROOT=/tmp/packstar-case-preflight-${HW_LABEL}-${SLURM_JOB_ID}
HOST_BUDGET_BYTES=${HOST_BUDGET_BYTES:-450971566080}
CASE_SLA_HOURS=${CASE_SLA_HOURS:-336}
SOFT_STATE_HOURS=${SOFT_STATE_HOURS:-24}
SAFETY_FACTOR=${SAFETY_FACTOR:-1.5}
INITIAL_EXACT_MILLIS=${INITIAL_EXACT_MILLIS:-120000}
FINAL_EXACT_MILLIS=${FINAL_EXACT_MILLIS:-300000}
FINAL_MAX_STATES=${FINAL_MAX_STATES:-3}

cd "$REPO"
mkdir -p "$RUN_ROOT/emat_cache"
cp -a "$CACHE_ROOT/3bua_fennario" "$RUN_ROOT/emat_cache/"
cp -a "$CACHE_ROOT/3k3q_fennario" "$RUN_ROOT/emat_cache/"

CP=$(<"$CP_FILE")
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx80g
    -Xms4g
    -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$RUN_ROOT/dp_mmap"
    -Dpackstar.decomp.strategy=adaptive
    -Dpackstar.decomp.weightedHicks.restarts=24
    -Dpackstar.decomp.weightedHicks.randomMoves=32
    -Dpackstar.decomp.exactImprove.maxDrop=1
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
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=0.5
    -Dosprey.bench.method=packstar
    -Dosprey.bench.packstarPreflightOnly=true
    -Dosprey.bench.numCPUs=16
    -Dosprey.bench.outputDir="$RUN_ROOT"
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

echo "preflight_job=$SLURM_JOB_ID hardware=$HW_LABEL gpuRatePerGpu=$GPU_RATE_PER_GPU gpuCount=$GPU_COUNT oocRate=$OOC_RATE gpuBudgetBytes=$GPU_BUDGET_BYTES safetyFactor=$SAFETY_FACTOR softStateHours=$SOFT_STATE_HOURS caseSlaHours=$CASE_SLA_HOURS start=$(date --iso-8601=seconds)"
echo "git_head=$(git rev-parse HEAD)"
sha256sum \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/BranchDpAdmission.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/BranchDpBackend.class \
    build/classes/java/main/edu/duke/cs/osprey/packstar/PackStarAdmissionDecision.class \
    build/classes/java/main/edu/duke/cs/osprey/packstar/PackStarCasePreflight.class \
    build/classes/java/test/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.class

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

echo "BEGIN case=3k3q hardware=$HW_LABEL time=$(date --iso-8601=seconds)"
if /usr/bin/time -f "PREFLIGHT_WALL hardware=$HW_LABEL design=3k3q elapsed=%e maxRssKiB=%M" \
    "$JAVA" "${JARGS[@]}" \
    -Dosprey.bench.designId=3k3q_fennario \
    -Dosprey.bench.pdbPath="$PDB_ROOT/3k3q/3k3q.min.reduce.renum.pdb" \
    '-Dosprey.bench.mutable=C444;C440' \
    '-Dosprey.bench.flexible=B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439' \
    -cp "$CP" "$MAIN"; then
    echo "END case=3k3q hardware=$HW_LABEL status=admitted time=$(date --iso-8601=seconds)"
else
    rc=$?
    status=1
    echo "END case=3k3q hardware=$HW_LABEL status=rejected-or-failed rc=$rc time=$(date --iso-8601=seconds)"
fi

echo "preflight_job=$SLURM_JOB_ID hardware=$HW_LABEL end=$(date --iso-8601=seconds) status=$status"
exit "$status"
