#!/bin/bash
#SBATCH --job-name=ra_scout
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_scout_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/ra_scout_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Region-Atom BranchMARK* SCOUT — diagnostic only.
# - Phase 0 emat-only pre-min L_R(b)/U_R(b) tables for each region atom (R, B=N_G(R)\R)
# - Phase 1 queue scout (record leaves before any minimization)
# - Phase 2 what-if epsilon-drop projections for delta = whatIfDeltas
# Does NOT certify a pfunc — exits early when SCOUT_MAX_LEAVES leaves observed.
# Wired to BranchMARK* SPARSE pfunc; no GPU required (Phase 0 is emat-only).

set -u

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
echo "  Start:      $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Isolate gradle caches per-job in case multiple scouts run in parallel.
export GRADLE_USER_HOME="/tmp/lz280-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE="/tmp/lz280-projcache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE"

echo "Compiling (isolated)..."
./gradlew --project-cache-dir "$PROJECT_CACHE" compileJava compileTestJava --no-daemon 2>&1 | tail -5
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "COMPILATION FAILED"
    rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
    exit 1
fi
echo "Compilation OK."
echo ""

./gradlew --project-cache-dir "$PROJECT_CACHE" test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testBranchMARKStarSparseOnly" \
    --no-daemon --rerun-tasks \
    -Dorg.gradle.jvmargs="-Xmx48g" \
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
