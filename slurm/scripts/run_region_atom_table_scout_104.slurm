#!/bin/bash
#SBATCH --job-name=ra_scout104
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_scout104_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_scout104_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-[01-06]
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Region-Atom BranchMARK* SCOUT — diagnostic only.
# - Phase 0 emat-only pre-min L_R(b)/U_R(b) tables for each region atom (R, B=N_G(R)\R)
# - Phase 1 queue scout (record leaves before any minimization)
# - Phase 2 what-if epsilon-drop projections for delta = whatIfDeltas
# Does NOT certify a pfunc — exits early when SCOUT_MAX_LEAVES leaves observed.
# Wired to BranchMARK* SPARSE pfunc; no GPU required (Phase 0 is emat-only).

set -u

# Constrain JVM common ForkJoinPool to match osprey thread count so that
# pending-edge lookahead (IntStream.parallel) uses all available CPUs.
export _JAVA_OPTIONS="-Djava.util.concurrent.ForkJoinPool.common.parallelism=104"

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
BUDGET=${BUDGET:-0.5}
ENERGY_MODE=${ENERGY_MODE:-sparse}
REGIONS=${REGIONS:-"7;6,7,8;0,2,6,7,8"}
SCOUT_MAX_LEAVES=${SCOUT_MAX_LEAVES:-2000}
TOP_CELLS=${TOP_CELLS:-12}
TABLE_MAX_JOBS=${TABLE_MAX_JOBS:-50000000}
WHATIF_DELTAS=${WHATIF_DELTAS:-"0.1,0.3,0.5"}
SCOUT_ONLY=${SCOUT_ONLY:-true}
# Phase 3 certify: ";"-separated atom names in canonical form (sorted comma-list),
# e.g. "7;6,7,8". Empty by default (scout only, no CCD oracle).
CERTIFY=${CERTIFY:-""}
# Phase 4 DP factor replacement (per-node bound tightening using the certified table).
# Only has any effect when CERTIFY is non-empty.
CERTIFY_USE_DP=${CERTIFY_USE_DP:-false}
# Top-K boundary cells (by Phase 0 emat gap) to certify. 0 = all cells.
# For atoms with many boundary cells, set this to ~100-1000 to keep Phase 3 feasible.
CERTIFY_TOPK=${CERTIFY_TOPK:-0}
# Edge selection strategy: lambda_states (default) | contraction | contraction_per_state
EDGE_SELECTION=${EDGE_SELECTION:-lambda_states}
# Root split strategy: legacy (default, fixed split) | lookahead (lookahead-driven, exposes multi-pending edges)
ROOT_SPLIT=${ROOT_SPLIT:-legacy}
ROOT_SPLIT_MAX_FSET=${ROOT_SPLIT_MAX_FSET:-2}

echo "=============================================="
echo "Region-Atom BranchMARK* SCOUT"
echo "  Job:        $SLURM_JOB_ID"
echo "  Node:       $SLURM_NODELIST"
echo "  CPUs:       $SLURM_CPUS_PER_TASK"
echo "  NUM_FLEX:   $NUM_FLEXIBLE"
echo "  EPSILON:    $EPSILON"
echo "  BUDGET:     $BUDGET (residual budget for sparse graph)"
echo "  ENERGY:     $ENERGY_MODE"
echo "  REGIONS:    $REGIONS"
echo "  MAX_LEAVES: $SCOUT_MAX_LEAVES"
echo "  TOP_CELLS:  $TOP_CELLS"
echo "  MAX_JOBS:   $TABLE_MAX_JOBS"
echo "  DELTAS:     $WHATIF_DELTAS"
echo "  SCOUT_ONLY: $SCOUT_ONLY"
echo "  CERTIFY:    $CERTIFY"
echo "  USE_DP:     $CERTIFY_USE_DP"
echo "  TOPK:       $CERTIFY_TOPK"
echo "  EDGE_SEL:   $EDGE_SELECTION"
echo "  ROOT_SPLIT: $ROOT_SPLIT (maxFset=$ROOT_SPLIT_MAX_FSET)"
echo "  Start:      $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean stale NFS lock files from prior gradle runs ('.nfsXXXX' silly-renames left
# behind when a process held an open file handle when files were deleted).
# Without this, gradle :processResources / :test fails with "Unable to delete directory".
echo "Cleaning stale NFS lock files..."
find /home/users/lz280/IdeaProjects/OSPREY3/build -name ".nfs*" -delete 2>/dev/null
# Some lock dirs can't be partially deleted; remove the offending subtrees outright.
rm -rf /home/users/lz280/IdeaProjects/OSPREY3/build/resources/main/config 2>/dev/null
rm -rf /home/users/lz280/IdeaProjects/OSPREY3/build/test-results/test/binary 2>/dev/null

# Isolate gradle caches per-job in case multiple scouts run in parallel.
export GRADLE_USER_HOME="/tmp/lz280-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE="/tmp/lz280-projcache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE"

echo "Compiling (isolated)..."
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
echo "Region-atom scout summary (parseable prefixes):"
echo "  [REGION_ATOM]            init / early-exit notes"
echo "  [REGION_ATOM_SPEC]       atom specs (R, B, table dims, totalJobs)"
echo "  [REGION_ATOM_TABLE]      emat-only Phase 0 table totals"
echo "  [REGION_ATOM_BOUNDARY]   top-K boundary cells (Phase 0 + queue)"
echo "  [REGION_ATOM_LOCAL]      top-K per-(b, x_R) local cells"
echo "  [REGION_ATOM_SCOUT]      Phase 1 queue scout summary"
echo "  [REGION_ATOM_WHATIF]     Phase 2 epsilon-drop simulation"
echo "  [REGION_ATOM_CERT]       Phase 3 certified-table aggregates"
echo "  [REGION_ATOM_CERT_BOUNDARY]  top-K certified boundary cells"
echo "  [REGION_ATOM_DP]         Phase 4 DP factor replacement counters"
echo ""
echo "End: $(date)"
rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
