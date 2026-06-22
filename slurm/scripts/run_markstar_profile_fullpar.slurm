#!/bin/bash
#SBATCH --job-name=msprof_fp
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/markstar_profile_fullpar_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/markstar_profile_fullpar_%j.err
#SBATCH --time=12:00:00
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
LEAF_PROFILE_DETAIL=${LEAF_PROFILE_DETAIL:-100000}
LEAF_PROFILE_INTERVAL=${LEAF_PROFILE_INTERVAL:-1}
NUM_CPUS=${NUM_CPUS:-104}

echo "=============================================="
echo "MARK* ONLY profile (fullParallelFromStart=true)"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "  NUM_FLEXIBLE=$NUM_FLEXIBLE"
echo "  EPSILON=$EPSILON"
echo "  LEAF_PROFILE_DETAIL=$LEAF_PROFILE_DETAIL"
echo "  LEAF_PROFILE_INTERVAL=$LEAF_PROFILE_INTERVAL"
echo "  LEAF_PROFILE_ALL_CONFS=true"
echo "  FULL_PARALLEL_FROM_START=true"
echo "  NUM_CPUS=$NUM_CPUS"
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
       -Xms8g -Xmx400g"

$JAVA $JARGS \
    -Dbranchdp.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchdp.test.epsilon="$EPSILON" \
    -Dmarkstar.leafProfile=true \
    -Dmarkstar.leafProfile.allConfs=true \
    -Dmarkstar.leafProfile.detail="$LEAF_PROFILE_DETAIL" \
    -Dmarkstar.leafProfile.interval="$LEAF_PROFILE_INTERVAL" \
    -Dmarkstar.fullParallelFromStart=true \
    -Dosprey.branchdp.numCpus="$NUM_CPUS" \
    -cp "$CP" \
    edu.duke.cs.osprey.markstar.TestBranchMARKStar markstar

echo ""
echo "End: $(date)"
