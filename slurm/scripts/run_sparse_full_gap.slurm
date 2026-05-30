#!/bin/bash
#SBATCH --job-name=sparse_full_gap
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/sparse_full_gap_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/sparse_full_gap_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=compsci
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "Sparse vs Full Partition Function Gap"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

NUM_FLEXIBLE=${NUM_FLEXIBLE:-3}
MAX_CONFS=${MAX_CONFS:-1000000}

echo ""
echo "Benchmark settings:"
echo "  NUM_FLEXIBLE=$NUM_FLEXIBLE"
echo "  MAX_CONFS=$MAX_CONFS"

echo ""
echo "Compiling tests..."
./gradlew compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running compareSparseFullPartitionFunctionGap"
echo "=========================================="
./gradlew test \
    --tests "edu.duke.cs.osprey.markstar.bench.SparseFullCheck.compareSparseFullPartitionFunctionGap" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx28g" \
    -Dbranchmarkstar.sparseFull.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchmarkstar.sparseFull.maxConfs="$MAX_CONFS" \
    2>&1

echo ""
echo "=========================================="
echo "Key results (from XML report):"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.markstar.bench.SparseFullCheck.xml
if [ -f "$XML" ]; then
    grep -E 'SPARSE/FULL|log10Z|K\\* log10|sparse-full|InteractionGraph dual cutoff' "$XML"
else
    echo "No report found"
fi

echo ""
echo "End time: $(date)"
