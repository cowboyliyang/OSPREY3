#!/bin/bash
#SBATCH --job-name=edge_sel_cmp
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/edge_selection_cmp_grisman_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/edge_selection_cmp_grisman_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-04
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "BranchMARK* Edge Selection Comparison (grisman)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
LOOKAHEAD_MAX_STATES=${LOOKAHEAD_MAX_STATES:-2000}
PARALLEL_LOOKAHEAD=${PARALLEL_LOOKAHEAD:-true}
PARALLEL_INTERNAL=${PARALLEL_INTERNAL:-false}
PARALLEL_ENUMERATION=${PARALLEL_ENUMERATION:-false}
EDGE_SELECTION=${EDGE_SELECTION:-contractionPerState}
TEST_METHOD=${TEST_METHOD:-testFlatSumVsOriginalMARKStar}
LEAF_PROFILE_DETAIL=${LEAF_PROFILE_DETAIL:-20}
LEAF_PROFILE_INTERVAL=${LEAF_PROFILE_INTERVAL:-100}

echo ""
echo "Benchmark settings:"
echo "  NUM_FLEXIBLE=$NUM_FLEXIBLE"
echo "  EPSILON=$EPSILON"
echo "  LOOKAHEAD_MAX_STATES=$LOOKAHEAD_MAX_STATES"
echo "  PARALLEL_LOOKAHEAD=$PARALLEL_LOOKAHEAD"
echo "  PARALLEL_INTERNAL=$PARALLEL_INTERNAL (DEPRECATED, forced off)"
echo "  PARALLEL_ENUMERATION=$PARALLEL_ENUMERATION (DEPRECATED, forced off)"
echo "  EDGE_SELECTION=$EDGE_SELECTION"
echo "  LEAF_PROFILE_DETAIL=$LEAF_PROFILE_DETAIL"
echo "  LEAF_PROFILE_INTERVAL=$LEAF_PROFILE_INTERVAL"
echo "  LEAF_PROFILE_ALL_CONFS=true"

echo ""
echo "Cleaning Gradle locks and daemons..."
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1 | tail -30
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation successful."

echo ""
echo "=========================================="
echo "Running testFlatSumVsOriginalMARKStar"
echo "BranchMARK* edge selection fixed to $EDGE_SELECTION"
echo "=========================================="
echo ""
TEST_JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Dbranchmarkstar.test.numFlexible=$NUM_FLEXIBLE -Dbranchmarkstar.test.epsilon=$EPSILON -Dbranchmarkstar.edgeSelection=$EDGE_SELECTION -Dbranchmarkstar.edgeSelection.maxExactStates=$LOOKAHEAD_MAX_STATES -Dbranchmarkstar.edgeSelection.parallelLookahead=$PARALLEL_LOOKAHEAD -Dbranchmarkstar.parallel.internal=$PARALLEL_INTERNAL -Dbranchmarkstar.parallel.enumeration=$PARALLEL_ENUMERATION -Dmarkstar.leafProfile=true -Dmarkstar.leafProfile.allConfs=true -Dmarkstar.leafProfile.detail=$LEAF_PROFILE_DETAIL -Dmarkstar.leafProfile.interval=$LEAF_PROFILE_INTERVAL"
echo "  TEST_JAVA_TOOL_OPTIONS=$TEST_JAVA_TOOL_OPTIONS"
echo ""
JAVA_TOOL_OPTIONS="$TEST_JAVA_TOOL_OPTIONS" ./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.$TEST_METHOD" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx48g" \
    -Dbranchmarkstar.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchmarkstar.test.epsilon="$EPSILON" \
    -Dbranchmarkstar.edgeSelection="$EDGE_SELECTION" \
    -Dbranchmarkstar.edgeSelection.maxExactStates="$LOOKAHEAD_MAX_STATES" \
    -Dbranchmarkstar.edgeSelection.parallelLookahead="$PARALLEL_LOOKAHEAD" \
    -Dbranchmarkstar.parallel.internal="$PARALLEL_INTERNAL" \
    -Dbranchmarkstar.parallel.enumeration="$PARALLEL_ENUMERATION" \
    -Dmarkstar.leafProfile=true \
    -Dmarkstar.leafProfile.allConfs=true \
    -Dmarkstar.leafProfile.detail="$LEAF_PROFILE_DETAIL" \
    -Dmarkstar.leafProfile.interval="$LEAF_PROFILE_INTERVAL" \
    2>&1

echo ""
echo "=========================================="
echo "Key results (from XML report):"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestBranchMARKStar.xml
if [ -f "$XML" ]; then
    grep -E 'LEAF_PROFILE|LEAF_PROFILE_SUMMARY|COMPARISON|Original MARK\*|BranchMARK\*|Speedup|BranchMARK\*:|epsilon=|enum steps|drill-downs|parallelCalls' "$XML" | head -300
else
    echo "No report found"
fi

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
