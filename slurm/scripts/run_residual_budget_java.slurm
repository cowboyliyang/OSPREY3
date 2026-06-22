#!/bin/bash
#SBATCH --job-name=resB
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_bb%x_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_bb%x_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --exclusive
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -u

BUDGET=${BUDGET:?must set BUDGET}
NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
ENERGY_MODE=${ENERGY_MODE:-sparse}
USE_CORRECTIONS=${USE_CORRECTIONS:-false}

echo "=============================================="
echo "BranchMARK* residualBudget=$BUDGET (eps=$EPSILON, numFlex=$NUM_FLEXIBLE)"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
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
       -Xms4g -Xmx72g"

$JAVA $JARGS \
    -Dbranchdp.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchdp.test.epsilon="$EPSILON" \
    -Dbranchdp.cutoff.strategy=residual_budget \
    -Dbranchdp.cutoff.residualBudget="$BUDGET" \
    -Dbranchdp.cutoff.keepConnected=true \
    -Dbranchdp.energyMode="$ENERGY_MODE" \
    -Dbranchdp.useHigherOrderCorrections="$USE_CORRECTIONS" \
    -cp "$CP" \
    edu.duke.cs.osprey.markstar.TestBranchMARKStar branchsparse

echo ""
echo "End: $(date)"
