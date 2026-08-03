#!/bin/bash
# 3bua with N all-20 MUTABLE positions (no flexible) -> preflight predictedCaseHours.
# Conf space differs from cached design, so the emat is rebuilt (CPU).
# Run: env N=6 BUDGET=0.5, or job array 0..11 (N in {6,7,8} x budget {0.5,1,2,3}).
#SBATCH --job-name=3bua_frontier
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err

set -euo pipefail
MUTLIST=(C447 C450 D638 D641 D642 D645 D623 D619 D616 D620 C417 C414 C410 C421 C452 C443 C444 C446)  # 3bua full 18-pos interface cluster, tightest-first
RB=${RB:-2}
if [ -n "${N:-}" ]; then
    n="$N"
else
    n=$(( ${SLURM_ARRAY_TASK_ID:?need N or SLURM_ARRAY_TASK_ID} + 2 ))   # array 0..16 -> N=2..18
fi
budget="$RB"
mutable=$(IFS=';'; echo "${MUTLIST[*]:0:$n}")   # pure all-20, flexible stays empty

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB=/usr/xtmp/lz280/dance_bench/pdbs_prepped/3bua/3bua.min.reduce.renum.pdb
RUN_ROOT=/usr/xtmp/lz280/packstar_production/all20/3bua_frontier_N${n}_rb${budget}_${SLURM_JOB_ID}
POLICY_DIR=/usr/xtmp/lz280/packstar_production/preflight_policies
POLICY_OUT=$POLICY_DIR/3bua_frontier_N${n}_rb${budget}.tsv
THREADS=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$RUN_ROOT/emat_cache/3bua" "$POLICY_DIR"

cd "$REPO"; CP=$(<"$CP_FILE")
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx100g -Xms4g -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$RUN_ROOT/dp_mmap"
    -Dpackstar.decomp.strategy=adaptive
    -Dpackstar.decomp.weightedHicks.restarts=32
    -Dpackstar.decomp.weightedHicks.randomMoves=32
    -Dpackstar.decomp.exactImprove.maxDrop=2
    -Dpackstar.decomp.exactImprove.maxMillis=120000
    -Dpackstar.rootSplit=predicted
    -Dpackstar.rootSplit.hostBudgetBytes=450971566080
    -Dpackstar.rootSplit.gpuBudgetBytes=23085449216
    -Dpackstar.dp.gpu=true
    -Dpackstar.dp.gpu.childSliceMaxBytes=2147483648
    -Dpackstar.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dpackstar.admission.gpuWorkPerSecondPerGpu=614000000
    -Dpackstar.admission.gpuCount=8
    -Dpackstar.admission.oocBytesPerSecond=170000000
    -Dpackstar.admission.safetyFactor=1.5
    -Dpackstar.admission.softStateHours=24
    -Dpackstar.admission.caseSlaHours=336
    -Dpackstar.admission.finalMaxStates=32
    -Dpackstar.admission.finalExactMaxMillis=120000
    -Dpackstar.admission.maxRounds=3
    -Dpackstar.admission.previewThreads="$THREADS"
    -Dpackstar.admission.policyOut="$POLICY_OUT"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget="$budget"
    -Dosprey.bench.method=packstar
    -Dosprey.bench.packstarPreflightOnly=true
    -Dosprey.bench.numCPUs="$THREADS"
    -Dosprey.bench.outputDir="$RUN_ROOT"
    -Dosprey.bench.designId=3bua
    -Dosprey.bench.pdbPath="$PDB"
    "-Dosprey.bench.mutable=$mutable"
    -Dosprey.bench.flexible=
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
echo "3bua_frontier N=$n rb=$budget mutable=$mutable threads=$THREADS node=$SLURMD_NODENAME start=$(date --iso-8601=seconds)"
status=0
if /usr/bin/time -f "ALL20_WALL N=$n budget=$budget elapsed=%e maxRssKiB=%M" \
    "$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"; then
    echo "END N=$n budget=$budget status=admitted time=$(date --iso-8601=seconds)"
else
    rc=$?; status=1; echo "END N=$n budget=$budget status=rejected-or-failed rc=$rc time=$(date --iso-8601=seconds)"
fi
rm -rf "$RUN_ROOT" 2>/dev/null || true
exit "$status"
