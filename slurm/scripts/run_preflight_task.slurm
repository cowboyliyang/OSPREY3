#!/bin/bash
# One PACK* admission preflight for (design, residualBudget), reusing the
# existing MARK*-format emats (export.* copied to packstar.* per task).
# Runs as a plain job (DESIGN=.. BUDGET=.. env) or a job array (0..113).
#SBATCH --job-name=pf_task
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err

set -euo pipefail
CSV=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv
BUDGETS=(0.5 1 2)

# ---- resolve design + budget (env override, else array index) ----
if [ -n "${DESIGN:-}" ] && [ -n "${BUDGET:-}" ]; then
    row=$(awk -F, -v d="$DESIGN" 'NR>2 && $1==d {print; exit}' "$CSV")
    budget="$BUDGET"
else
    tid=${SLURM_ARRAY_TASK_ID:?need SLURM_ARRAY_TASK_ID or DESIGN+BUDGET}
    didx=$((tid/3)); bidx=$((tid%3)); budget=${BUDGETS[$bidx]}
    row=$(awk -F, -v n=$((didx+1)) 'NR>2 {c++; if(c==n){print; exit}}' "$CSV")
fi
design_id=$(printf '%s' "$row" | awk -F, '{print $1}')
pdb=$(printf '%s' "$row" | awk -F, '{print $2}')
mutable=$(printf '%s' "$row" | awk -F, '{print $6}')
flexible=$(printf '%s' "$row" | awk -F, '{print $7}')
[ -n "$design_id" ] || { echo "no design resolved"; exit 2; }

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
SRC=/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id
RUN_ROOT=/usr/xtmp/lz280/packstar_production/preflight/${design_id}_rb${budget}_${SLURM_JOB_ID}
EDIR=$RUN_ROOT/emat_cache/$design_id
POLICY_DIR=/usr/xtmp/lz280/packstar_production/preflight_policies
POLICY_OUT=$POLICY_DIR/${design_id}_rb${budget}.tsv
THREADS=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$EDIR" "$POLICY_DIR"

# ---- reuse existing emats: copy export.* -> packstar.* (isolated per task) ----
for st in complex ligand protein; do
    cp "$SRC/export.$st.rigid.dat" "$EDIR/packstar.$st.rigid.dat"
    cp "$SRC/export.$st.min.dat"   "$EDIR/packstar.$st.minimizing.dat"
done

cd "$REPO"
CP=$(<"$CP_FILE")
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx80g -Xms4g -XX:-UseSuperWord
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
    -Dpackstar.admission.finalMaxStates=64
    -Dpackstar.admission.finalExactMaxMillis=180000
    -Dpackstar.admission.maxRounds=4
    -Dpackstar.admission.previewThreads="$THREADS"
    -Dpackstar.admission.policyOut="$POLICY_OUT"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget="$budget"
    -Dosprey.bench.method=packstar
    -Dosprey.bench.packstarPreflightOnly=true
    -Dosprey.bench.numCPUs="$THREADS"
    -Dosprey.bench.outputDir="$RUN_ROOT"
    -Dosprey.bench.designId="$design_id"
    -Dosprey.bench.pdbPath="$PDB_ROOT/$pdb/$pdb.min.reduce.renum.pdb"
    "-Dosprey.bench.mutable=$mutable"
    "-Dosprey.bench.flexible=$flexible"
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

echo "pf_task design=$design_id pdb=$pdb budget=$budget threads=$THREADS policyOut=$POLICY_OUT node=$SLURMD_NODENAME start=$(date --iso-8601=seconds)"
echo "mutable=$mutable"; echo "flexible=$flexible"
status=0
if /usr/bin/time -f "PREFLIGHT_WALL design=$design_id budget=$budget elapsed=%e maxRssKiB=%M" \
    "$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"; then
    echo "END design=$design_id budget=$budget status=admitted time=$(date --iso-8601=seconds)"
else
    rc=$?; status=1
    echo "END design=$design_id budget=$budget status=rejected-or-failed rc=$rc time=$(date --iso-8601=seconds)"
fi
rm -rf "$RUN_ROOT" 2>/dev/null || true   # policy dump is kept in POLICY_DIR
exit "$status"
