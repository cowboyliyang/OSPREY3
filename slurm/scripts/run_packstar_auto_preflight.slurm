#!/bin/bash
# Usage: sbatch --account=grisman --export=ALL,DESIGN=1a0r this-file
#SBATCH --job-name=ps-auto-preflight
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${DESIGN:?Supply a design ID from design_specs_prepped.csv}"
AUTO_PREFLIGHT_CLASSPATH=${AUTO_PREFLIGHT_CLASSPATH:-/usr/xtmp/lz280/packstar_auto_preflight/current/test_classpath.txt}
test -s "$AUTO_PREFLIGHT_CLASSPATH"
CSV=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv
row=$(awk -F, -v d="$DESIGN" 'NR>2 && $1==d {print; exit}' "$CSV")
test -n "$row"
IFS=, read -r design_id pdb description num_seq num_pos mutable flexible remainder <<< "$row"
RUN_ROOT=/usr/xtmp/lz280/packstar_auto_preflight/${design_id}_${SLURM_JOB_ID}
mkdir -p "$RUN_ROOT/tmp" "$RUN_ROOT/emat_cache/$design_id"
export TMPDIR="$RUN_ROOT/tmp"
export CUDA_CACHE_PATH="$RUN_ROOT/cuda-cache"
THREADS=${SLURM_CPUS_PER_TASK:-8}
# Heap sizing follows the allocation; no target GPU model or throughput is supplied.
MEM_MIB=${SLURM_MEM_PER_NODE:-65536}
HEAP_MIB=$((MEM_MIB * 3 / 4))
CP=$(<"$AUTO_PREFLIGHT_CLASSPATH")
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
SOURCE_CACHE=/usr/xtmp/lz280/bench_comparison/results/emat_cache/$design_id
# Existing small per-design energy matrices only; absence triggers normal EMAT construction.
for state in complex protein ligand; do
    if [ -f "$SOURCE_CACHE/export.$state.rigid.dat" ] && [ -f "$SOURCE_CACHE/export.$state.min.dat" ]; then
        cp "$SOURCE_CACHE/export.$state.rigid.dat" "$RUN_ROOT/emat_cache/$design_id/packstar.$state.rigid.dat"
        cp "$SOURCE_CACHE/export.$state.min.dat" "$RUN_ROOT/emat_cache/$design_id/packstar.$state.minimizing.dat"
    fi
done
cd "$RUN_ROOT"
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    "-Xmx${HEAP_MIB}m" -XX:-UseSuperWord "-Djava.io.tmpdir=$TMPDIR"
    -Dpackstar.admission.mode=auto
    "-Dpackstar.admission.policyOut=$RUN_ROOT/policy.tsv"
    -Dpackstar.dp.cache=false -Dbranchdp.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap "-Dbranchdp.dp.mmap.dir=$RUN_ROOT/dp_mmap"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET "-Dbranchdp.cutoff.residualBudget=${BUDGET:-0.5}"
    -Dosprey.bench.method=packstar -Dosprey.bench.packstarPreflightOnly=true
    "-Dosprey.bench.numCPUs=$THREADS" "-Dosprey.bench.outputDir=$RUN_ROOT"
    "-Dosprey.bench.designId=$design_id"
    "-Dosprey.bench.pdbPath=/usr/xtmp/lz280/dance_bench/pdbs_prepped/$pdb/$pdb.min.reduce.renum.pdb"
    "-Dosprey.bench.mutable=$mutable" "-Dosprey.bench.flexible=$flexible"
)
if [ -n "${POLICY_IN:-}" ]; then JARGS+=("-Dpackstar.admission.policyIn=$POLICY_IN"); fi
"$JAVA" "${JARGS[@]}" -cp "$CP" edu.duke.cs.osprey.markstar.bench.GenericPDBBench
echo "AUTO_PREFLIGHT_RESULTS=$RUN_ROOT"
