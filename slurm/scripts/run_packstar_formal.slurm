#!/bin/bash
# Formal PACK* K* run for one design at residualBudget=RB (default 2), on a GPU
# node. Reuses the design's emat (export.* == packstar.*, verified identical).
# Produces the sequence-ranking CSV. Run: env DESIGN=..., or job array 0..37.
#SBATCH --job-name=pf_formal
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err

set -euo pipefail
CSV=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv
if [ -n "${DESIGN:-}" ]; then
    row=$(awk -F, -v d="$DESIGN" 'NR>2 && $1==d {print; exit}' "$CSV")
else
    tid=${SLURM_ARRAY_TASK_ID:?need DESIGN or SLURM_ARRAY_TASK_ID}
    row=$(awk -F, -v n=$((tid+1)) 'NR>2 {c++; if(c==n){print; exit}}' "$CSV")
fi
design_id=$(printf '%s' "$row" | awk -F, '{print $1}')
pdb=$(printf '%s' "$row" | awk -F, '{print $2}')
mutable=$(printf '%s' "$row" | awk -F, '{print $6}')
flexible=$(printf '%s' "$row" | awk -F, '{print $7}')
RB=${RB:-2}

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
SRC=/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id
OUT_DIR=/usr/xtmp/lz280/packstar_production/formal_rb${RB}/${design_id}_${SLURM_JOB_ID}
EDIR=$OUT_DIR/emat_cache/$design_id
RANKDIR=/usr/xtmp/lz280/packstar_production/formal_rb${RB}_rankings
THREADS=${SLURM_CPUS_PER_TASK:-32}
mkdir -p "$EDIR" "$RANKDIR"
for st in complex ligand protein; do
    cp "$SRC/export.$st.rigid.dat" "$EDIR/packstar.$st.rigid.dat"
    cp "$SRC/export.$st.min.dat"   "$EDIR/packstar.$st.minimizing.dat"
done

cd "$REPO"; CP=$(<"$CP_FILE")
command -v module >/dev/null 2>&1 && { module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true; }
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader -i 0 2>/dev/null || true
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx280g -Xms4g -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.dp.mmap.dir="$OUT_DIR/dp_mmap"
    -Dbranchdp.decomp.strategy=weighted_hicks
    -Dbranchdp.decomp.weightedHicks.restarts=24
    -Dbranchdp.decomp.weightedHicks.randomMoves=32
    -Dbranchdp.dp.gpu=true
    -Dbranchdp.dp.gpu.childSliceMaxBytes=2147483648
    -Dbranchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dbranchdp.rootSplit=gpubytes
    -Dbranchdp.rootSplit.hostBudgetBytes=280000000000
    -Dbranchdp.rootSplit.gpuBudgetBytes=11000000000
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget="$RB"
    -Dosprey.bench.method=packstar
    -Dosprey.bench.numCPUs="$THREADS"
    -Dosprey.bench.outputDir="$OUT_DIR"
    -Dosprey.bench.designId="$design_id"
    -Dosprey.bench.pdbPath="$PDB_ROOT/$pdb/$pdb.min.reduce.renum.pdb"
    "-Dosprey.bench.mutable=$mutable"
    "-Dosprey.bench.flexible=$flexible"
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
echo "FORMAL design=$design_id rb=$RB node=$SLURMD_NODENAME gpus=${CUDA_VISIBLE_DEVICES:-?} start=$(date --iso-8601=seconds)"
echo "mutable=$mutable"; echo "flexible=$flexible"
status=0
if /usr/bin/time -f "FORMAL_WALL design=$design_id rb=$RB elapsed=%e maxRssKiB=%M" \
    "$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"; then
    csv="$OUT_DIR/${design_id}_packstar.csv"
    [ -f "$csv" ] && cp "$csv" "$RANKDIR/${design_id}.csv" && echo "ranking -> $RANKDIR/${design_id}.csv"
    echo "END design=$design_id rb=$RB status=ok time=$(date --iso-8601=seconds)"
else
    rc=$?; status=1; echo "END design=$design_id rb=$RB status=FAILED rc=$rc time=$(date --iso-8601=seconds)"
fi
rm -rf "$OUT_DIR/dp_mmap" 2>/dev/null || true
exit "$status"
