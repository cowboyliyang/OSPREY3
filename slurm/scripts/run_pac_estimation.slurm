#!/bin/bash
#SBATCH --job-name=pac_est
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/pac_estimation_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/pac_estimation_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=512G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --exclusive
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -u

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
PAC_SAMPLES=${PAC_SAMPLES:-500}
PAC_CONFIDENCE=${PAC_CONFIDENCE:-0.05}
PAC_RESIDUAL_BOUND=${PAC_RESIDUAL_BOUND:-0.5}
NUM_CPUS=${NUM_CPUS:-104}
DP_CACHE=${DP_CACHE:-false}
EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord"

echo "=============================================="
echo "PAC Partition Function Estimation"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "  NUM_FLEXIBLE=$NUM_FLEXIBLE"
echo "  EPSILON=$EPSILON"
echo "  PAC_SAMPLES=$PAC_SAMPLES"
echo "  PAC_CONFIDENCE=$PAC_CONFIDENCE  (=> $(echo "scale=0; (1-$PAC_CONFIDENCE)*100" | bc)% confidence)"
echo "  PAC_RESIDUAL_BOUND=$PAC_RESIDUAL_BOUND kcal/mol"
echo "  NUM_CPUS=$NUM_CPUS"
echo "  DP_CACHE=$DP_CACHE"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java

# Build classpath from already-compiled classes + dependency jars (no gradle)
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do
    CP="$CP:$jar"
done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do
    CP="$CP:$jar"
done
echo "Classpath entries: $(echo $CP | tr ':' '\n' | wc -l)"
echo ""

JARGS="--add-opens java.base/java.util=ALL-UNNAMED \
       --add-opens java.base/java.lang=ALL-UNNAMED \
       --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
       -Xms8g -Xmx400g \
       $EXTRA_JVM_ARGS"

$JAVA $JARGS \
    -Dbranchmarkstar.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchmarkstar.test.epsilon="$EPSILON" \
    -Dbranchmarkstar.usePAC=true \
    -Dbranchmarkstar.pac.residualBound="$PAC_RESIDUAL_BOUND" \
    -Dbranchmarkstar.dp.cache="$DP_CACHE" \
    -Dbranchmarkstar.pac.samples="$PAC_SAMPLES" \
    -Dbranchmarkstar.pac.confidence="$PAC_CONFIDENCE" \
    -Dosprey.branchmarkstar.numCpus="$NUM_CPUS" \
    -cp "$CP" \
    edu.duke.cs.osprey.markstar.TestBranchMARKStar branchsparse

echo ""
echo "End: $(date)"
