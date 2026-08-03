#!/bin/bash
#SBATCH --job-name=3k3q_fennario_full
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-03
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=950G
#SBATCH --gres=gpu:a5000:8
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/3k3q_fennario_full_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/3k3q_fennario_full_%j.err

set -euo pipefail

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
OUT_DIR=/usr/xtmp/lz280/packstar_production/3k3q_fennario_${SLURM_JOB_ID}
BUILD_ROOT=$REPO/slurm/artifacts/packstar_gpu_ooc_prod_20260710_r1

cd "$REPO"
mkdir -p "$OUT_DIR"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi
command -v module >/dev/null 2>&1 && { module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true; }

echo "node=$SLURM_NODELIST gpus=${CUDA_VISIBLE_DEVICES:-unset} start=$(date)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader

sha256sum "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuFullDP.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/BranchDpBackend.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/BranchDecomposition.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/DPTable.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/MappedDPTable.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/branchdp/RootedTreeEdge.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/kstar/KStar\$ConfSpaceInfo.class" \
    "$BUILD_ROOT/build/classes/java/main/edu/duke/cs/osprey/packstar/PackStarPartitionFunction.class" \
    "$BUILD_ROOT/build/resources/main/gpuKernels/cuda/dp.bin"
DEPS_CP=$(cut -d: -f5- "$CP_FILE")
CP="$BUILD_ROOT/build/classes/java/main:$BUILD_ROOT/build/classes/java/test:$BUILD_ROOT/build/resources/main:$BUILD_ROOT/build/resources/test:$DEPS_CP"
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx850g
    -Xms4g
    -XX:-UseSuperWord
    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode=auto_mmap
    -Dbranchdp.dp.mmap.thresholdBytes=137438953472
    -Dbranchdp.dp.mmap.dir="$OUT_DIR/dp_mmap"
    -Dbranchdp.dp.mmap.skipInitialFill=true
    -Dbranchdp.decomp.strategy=weighted_hicks
    -Dbranchdp.decomp.weightedHicks.restarts=24
    -Dbranchdp.decomp.weightedHicks.randomMoves=32
    -Dbranchdp.dp.gpu=true
    -Dbranchdp.dp.gpu.childSliceMaxBytes=2147483648
    -Dbranchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
    -Dbranchdp.rootSplit=gpubytes
    -Dbranchdp.rootSplit.hostBudgetBytes=450971566080
    -Dbranchdp.rootSplit.gpuBudgetBytes=23085449216
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=0.5
    -Dosprey.bench.method=packstar
    -Dosprey.bench.numCPUs=100
    -Dosprey.bench.outputDir="$OUT_DIR"
    -Dosprey.bench.designId=3k3q_fennario
    -Dosprey.bench.pdbPath="$PDB_ROOT/3k3q/3k3q.min.reduce.renum.pdb"
    '-Dosprey.bench.mutable=C444;C440'
    '-Dosprey.bench.flexible=B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439'
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

echo "BEGIN 3k3q_fennario $(date)"
/usr/bin/time -f 'PROD_WALL design=3k3q_fennario elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"
echo "END 3k3q_fennario $(date)"
