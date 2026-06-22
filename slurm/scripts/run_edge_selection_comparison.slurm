#!/bin/bash
#SBATCH --job-name=edge_selection_cmp
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/edge_selection_cmp_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/edge_selection_cmp_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "BranchMARK* Edge Selection Comparison"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Optional overrides:
#   sbatch --export=ALL,NUM_FLEXIBLE=12,EPSILON=0.68 slurm/scripts/run_edge_selection_comparison.slurm
NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
LOOKAHEAD_MAX_STATES=${LOOKAHEAD_MAX_STATES:-2000}
PARALLEL_LOOKAHEAD=${PARALLEL_LOOKAHEAD:-true}

echo ""
echo "Benchmark settings:"
echo "  NUM_FLEXIBLE=$NUM_FLEXIBLE"
echo "  EPSILON=$EPSILON"
echo "  LOOKAHEAD_MAX_STATES=$LOOKAHEAD_MAX_STATES"
echo "  PARALLEL_LOOKAHEAD=$PARALLEL_LOOKAHEAD"

echo ""
echo "Cleaning Gradle locks and daemons..."
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation successful."

echo ""
echo "=========================================="
echo "Running testBranchEdgeSelectionComparison"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testBranchEdgeSelectionComparison" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx28g" \
    -Dbranchdp.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchdp.test.epsilon="$EPSILON" \
    -Dbranchdp.edgeSelection.maxExactStates="$LOOKAHEAD_MAX_STATES" \
    -Dbranchdp.edgeSelection.parallelLookahead="$PARALLEL_LOOKAHEAD" \
    2>&1

echo ""
echo "=========================================="
echo "Key results (from XML report):"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestBranchMARKStar.xml
if [ -f "$XML" ]; then
    grep -oP 'time="[^"]*"' "$XML" | head -5
    grep -E 'EDGE SELECTION|lambdaStates.*runtime|contraction.*runtime|Speedup|BranchMARK\\*: Edge selection|BranchMARK\\*: Edge selection stats|BranchMARK\\*: Numeric audit|BranchMARK\\*: Finished|BranchMARK\\*: [0-9]+ minimizations|K\\*:|epsilon=|enum steps|drill-downs' "$XML"
else
    echo "No report found"
fi

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
