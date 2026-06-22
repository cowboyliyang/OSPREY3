#!/bin/bash
#SBATCH --job-name=par_ablation
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/par_ablation_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/par_ablation_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --exclusive
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-04
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "BranchMARK* Pending-Edge Lookahead Sweep (grisman)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
NUM_CPUS=${NUM_CPUS:-104}
ROOT_SPLIT=${ROOT_SPLIT:-lookahead}
ROOT_SPLIT_MAX_FSET=${ROOT_SPLIT_MAX_FSET:-2}
LOOKAHEAD_MAX_PENDING_EDGES=${LOOKAHEAD_MAX_PENDING_EDGES:-4}
PENDING_EDGE_SWEEP=${PENDING_EDGE_SWEEP:-2,3,4,5,6,8,12}
LOOKAHEAD_MAX_STATES=${LOOKAHEAD_MAX_STATES:-2000}
EDGE_SELECTIONS=${EDGE_SELECTIONS:-contraction,contractionPerState}
PARALLEL_LOOKAHEAD=${PARALLEL_LOOKAHEAD:-true}
PARALLEL_INTERNAL=false
PARALLEL_ENUMERATION=false

# Constrain JVM common ForkJoinPool to match osprey thread count
# pending-edge lookahead uses IntStream.parallel(), which defaults to ForkJoinPool.common
export _JAVA_OPTIONS="-Djava.util.concurrent.ForkJoinPool.common.parallelism=$NUM_CPUS"

# Use a per-job Gradle home so this job does not contend with login-node or
# other Slurm jobs for ~/.gradle/caches/journal-1.
export GRADLE_USER_HOME="${GRADLE_USER_HOME:-${SLURM_TMPDIR:-/tmp}/osprey3-gradle-${USER}-${SLURM_JOB_ID}}"
mkdir -p "$GRADLE_USER_HOME"
if [ -d "$HOME/.gradle/wrapper" ] && [ ! -d "$GRADLE_USER_HOME/wrapper" ]; then
    mkdir -p "$GRADLE_USER_HOME"
    cp -a "$HOME/.gradle/wrapper" "$GRADLE_USER_HOME/"
fi

echo ""
echo "Settings: NUM_FLEXIBLE=$NUM_FLEXIBLE  EPSILON=$EPSILON  NUM_CPUS=$NUM_CPUS"
echo "BranchMARK*: ROOT_SPLIT=$ROOT_SPLIT  ROOT_SPLIT_MAX_FSET=$ROOT_SPLIT_MAX_FSET"
echo "BranchMARK*: EDGE_SELECTIONS=$EDGE_SELECTIONS"
echo "BranchMARK*: PENDING_EDGE_SWEEP=$PENDING_EDGE_SWEEP  LOOKAHEAD_MAX_STATES=$LOOKAHEAD_MAX_STATES"
echo "BranchMARK*: parallelLookahead=$PARALLEL_LOOKAHEAD"
echo "BranchMARK*: parallelInternal=$PARALLEL_INTERNAL (DEPRECATED, fixed off)"
echo "BranchMARK*: parallelEnumeration=$PARALLEL_ENUMERATION (DEPRECATED, fixed off)"
echo "ForkJoinPool.common.parallelism=$NUM_CPUS"
echo "GRADLE_USER_HOME=$GRADLE_USER_HOME"
echo ""

./gradlew --stop 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1 | tail -10
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation OK"

echo ""
echo "=========================================="
echo "Running testPendingEdgeLookaheadSweep"
echo "=========================================="
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.bench.PendingEdgeSweep.testPendingEdgeLookaheadSweep" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx48g" \
    -Dosprey.branchdp.numFlexible="$NUM_FLEXIBLE" \
    -Dosprey.branchdp.epsilon="$EPSILON" \
    -Dosprey.branchdp.numCpus="$NUM_CPUS" \
    -Dosprey.branchdp.rootSplit="$ROOT_SPLIT" \
    -Dosprey.branchdp.rootSplit.maxFset="$ROOT_SPLIT_MAX_FSET" \
    -Dosprey.branchdp.pendingEdgeSweep="$PENDING_EDGE_SWEEP" \
    -Dosprey.branchdp.edgeSelections="$EDGE_SELECTIONS" \
    -Dosprey.branchdp.edgeSelection.maxPendingEdges="$LOOKAHEAD_MAX_PENDING_EDGES" \
    -Dosprey.branchdp.edgeSelection.maxExactStates="$LOOKAHEAD_MAX_STATES" \
    -Dosprey.branchdp.edgeSelection.parallelLookahead="$PARALLEL_LOOKAHEAD" \
    -Dosprey.branchdp.parallel.internal="$PARALLEL_INTERNAL" \
    -Dosprey.branchdp.parallel.enumeration="$PARALLEL_ENUMERATION" \
    2>&1

echo ""
echo "=========================================="
echo "Pending-edge sweep results from XML report:"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.bench.PendingEdgeSweep.xml
if [ -f "$XML" ]; then
    grep -E 'PENDING EDGE LOOKAHEAD SWEEP|pendingEdgeSweep|BranchMARK\*: Edge selection|BranchMARK\*: Root split edge|pendingHistogram|parallelCalls|DEPRECATED' "$XML" | head -200
else
    echo "No XML report"
fi

echo ""
echo "End time: $(date)"
