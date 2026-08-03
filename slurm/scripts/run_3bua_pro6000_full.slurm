#!/bin/bash
#SBATCH --job-name=3bua_pro6000_full
#SBATCH --partition=compsci-gpu
#SBATCH --nodelist=compsci-cluster-fitz-47
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=1000G
#SBATCH --gres=gpu:rtx_pro_6000:3
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/3bua_pro6000_full_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/3bua_pro6000_full_%j.err

set -euo pipefail

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
OUT_DIR=/usr/xtmp/lz280/packstar_production/3bua_pro6000_${SLURM_JOB_ID}
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
    -Xmx900g
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
    -Dbranchdp.rootSplit.gpuBudgetBytes=91268055040
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=0.5
    -Dosprey.bench.method=packstar
    -Dosprey.bench.numCPUs=120
    -Dosprey.bench.outputDir="$OUT_DIR"
    -Dosprey.bench.designId=3bua_pro6000
    -Dosprey.bench.pdbPath="$PDB_ROOT/3bua/3bua.min.reduce.renum.pdb"
    '-Dosprey.bench.mutable=C447;C450'
    '-Dosprey.bench.flexible=D638;D641;D642;D645;D623;D619;D616;D620;C417;C414;C410;C421;C452;C443;C444;C446'
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

echo "BEGIN 3bua_pro6000 $(date)"
/usr/bin/time -f 'PROD_WALL design=3bua_pro6000 elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" -cp "$CP" "$MAIN"
echo "END 3bua_pro6000 $(date)"
