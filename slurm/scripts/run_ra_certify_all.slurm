#!/bin/bash
#SBATCH --job-name=ra_certall
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_certall_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_certall_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --nodelist=fennario-[01-06]
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Region-Atom BranchMARK* — CERTIFY ALL boundary cells + full search.
# Key changes vs ra_scout104:
#   - CERTIFY_TOPK=0 → certify ALL boundary cells (not just top-200)
#   - SCOUT_ONLY=false → run full K* search after certification
#   - CERTIFY_USE_DP=true → apply certified table during search
#   - contraction_per_state + lookahead root split
#   - 104 CPUs + 1 A5000 on fennario

set -u

export _JAVA_OPTIONS="-Djava.util.concurrent.ForkJoinPool.common.parallelism=104"

NUM_FLEXIBLE=10
EPSILON=0.68
BUDGET=0.5
ENERGY_MODE=sparse
REGIONS="7"
SCOUT_MAX_LEAVES=2000
TOP_CELLS=12
TABLE_MAX_JOBS=50000000
WHATIF_DELTAS="0.1,0.3,0.5"
SCOUT_ONLY=false
CERTIFY="7"
CERTIFY_USE_DP=true
CERTIFY_TOPK=0          # 0 = certify ALL boundary cells
EDGE_SELECTION=contraction_per_state
ROOT_SPLIT=lookahead
ROOT_SPLIT_MAX_FSET=2

echo "=============================================="
echo "Region-Atom BranchMARK* CERTIFY-ALL + SEARCH"
echo "  Job:        $SLURM_JOB_ID"
echo "  Node:       $SLURM_NODELIST"
echo "  CPUs:       $SLURM_CPUS_PER_TASK"
echo "  GPU:        A5000"
echo "  NUM_FLEX:   $NUM_FLEXIBLE"
echo "  EPSILON:    $EPSILON"
echo "  BUDGET:     $BUDGET"
echo "  ENERGY:     $ENERGY_MODE"
echo "  REGIONS:    $REGIONS"
echo "  MAX_LEAVES: $SCOUT_MAX_LEAVES"
echo "  CERTIFY:    $CERTIFY"
echo "  USE_DP:     $CERTIFY_USE_DP"
echo "  TOPK:       $CERTIFY_TOPK (0=ALL)"
echo "  EDGE_SEL:   $EDGE_SELECTION"
echo "  ROOT_SPLIT: $ROOT_SPLIT (maxFset=$ROOT_SPLIT_MAX_FSET)"
echo "  Start:      $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "Cleaning stale NFS lock files..."
find /home/users/lz280/IdeaProjects/OSPREY3/build -name ".nfs*" -delete 2>/dev/null
rm -rf /home/users/lz280/IdeaProjects/OSPREY3/build/resources/main/config 2>/dev/null
rm -rf /home/users/lz280/IdeaProjects/OSPREY3/build/test-results/test/binary 2>/dev/null

export GRADLE_USER_HOME="/tmp/lz280-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE="/tmp/lz280-projcache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE"

echo "Compiling..."
./gradlew --project-cache-dir "$PROJECT_CACHE" compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
    exit 1
fi
echo "Compilation OK."
echo ""

./gradlew --project-cache-dir "$PROJECT_CACHE" test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testBranchMARKStarSparseOnly" \
    --no-daemon --rerun-tasks \
    -Dorg.gradle.jvmargs="-Xmx200g" \
    -Dosprey.branchmarkstar.numCpus=104 \
    -Dbranchmarkstar.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchmarkstar.test.epsilon="$EPSILON" \
    -Dbranchmarkstar.energyMode="$ENERGY_MODE" \
    -Dbranchmarkstar.cutoff.strategy=residual_budget \
    -Dbranchmarkstar.cutoff.residualBudget="$BUDGET" \
    -Dbranchmarkstar.cutoff.keepConnected=true \
    -Dbranchmarkstar.regionAtom.enabled=true \
    -Dbranchmarkstar.regionAtom.scoutOnly="$SCOUT_ONLY" \
    -Dbranchmarkstar.regionAtom.regions="$REGIONS" \
    -Dbranchmarkstar.regionAtom.scout.maxLeaves="$SCOUT_MAX_LEAVES" \
    -Dbranchmarkstar.regionAtom.scout.topCells="$TOP_CELLS" \
    -Dbranchmarkstar.regionAtom.table.maxJobs="$TABLE_MAX_JOBS" \
    -Dbranchmarkstar.regionAtom.whatIfDeltas="$WHATIF_DELTAS" \
    -Dbranchmarkstar.regionAtom.certify="$CERTIFY" \
    -Dbranchmarkstar.regionAtom.certify.useDP="$CERTIFY_USE_DP" \
    -Dbranchmarkstar.regionAtom.certify.topK="$CERTIFY_TOPK" \
    -Dbranchmarkstar.edgeSelection="$EDGE_SELECTION" \
    -Dbranchmarkstar.rootSplit="$ROOT_SPLIT" \
    -Dbranchmarkstar.rootSplit.maxFset="$ROOT_SPLIT_MAX_FSET" \
    2>&1

echo ""
echo "Key metrics to compare vs ra_scout104_11791940 (topK=200):"
echo "  [REGION_ATOM_CERT]       certified table size + gap shrink"
echo "  [REGION_ATOM_DP]         applied vs missingBoundary"
echo "  BranchMARK*: N minimizations, epsilon=X  convergence speed"
echo ""
echo "End: $(date)"
rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
