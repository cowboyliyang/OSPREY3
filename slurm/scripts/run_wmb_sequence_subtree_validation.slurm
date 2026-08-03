#!/bin/bash
#SBATCH --job-name=wmbis_seq
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/wmbis_seq_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/wmbis_seq_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

STATE=${STATE:-complex}
IBOUNDS=${IBOUNDS:-1,2}
SAMPLES=${SAMPLES:-16,64}
SEED=${SEED:-1}
CONFIDENCE=${CONFIDENCE:-0.05}
LOCAL_CAP_MAX_ASSIGNMENTS=${LOCAL_CAP_MAX_ASSIGNMENTS:-5000000}
PROB_TAIL_QUANTILE=${PROB_TAIL_QUANTILE:-0.99}
MUTABLE_RESIDUES=${MUTABLE_RESIDUES:-A156,A164}
PROTEIN_WT_FLEX=${PROTEIN_WT_FLEX:-A172,A192,A193}
LIGAND_WT_FLEX=${LIGAND_WT_FLEX:-G649,G650,G651}
ASSIGNMENTS=${ASSIGNMENTS:-}

echo "=============================================="
echo "WMB-IS sequence-subtree validation"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "Java: $(java -version 2>&1 | head -n 1)"
echo "STATE=$STATE IBOUNDS=$IBOUNDS SAMPLES=$SAMPLES SEED=$SEED"
echo "LOCAL_CAP_MAX_ASSIGNMENTS=$LOCAL_CAP_MAX_ASSIGNMENTS PROB_TAIL_QUANTILE=$PROB_TAIL_QUANTILE"
echo "MUTABLE_RESIDUES=$MUTABLE_RESIDUES"
echo "PROTEIN_WT_FLEX=$PROTEIN_WT_FLEX"
echo "LIGAND_WT_FLEX=$LIGAND_WT_FLEX"
echo "ASSIGNMENTS=${ASSIGNMENTS:-root/unassigned}"
echo "=============================================="

CP="build/classes/java/main:build/classes/java/test:build/classes/kotlin/main:build/classes/kotlin/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do
    CP="$CP:$jar"
done
while IFS= read -r jar; do
    CP="$CP:$jar"
done < <(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null \
    | grep -v "onnxruntime" \
    | grep -v -- "-sources.jar" \
    | sort -u)

echo "Classpath entries: $(echo "$CP" | tr ':' '\n' | wc -l)"

echo "Compiling WMB and SPACK* validation sources..."
javac -cp "$CP" -d build/classes/java/main \
    src/main/java/edu/duke/cs/osprey/wmb/WmbModel.java \
    src/main/java/edu/duke/cs/osprey/wmb/WeightedMiniBucket.java \
    src/main/java/edu/duke/cs/osprey/wmb/MiniBucketBound.java
javac -cp "$CP" -d build/classes/java/main \
    $(find src/main/java/edu/duke/cs/osprey/spackstar -name '*.java')
javac -cp "$CP" -d build/classes/java/test \
    src/test/java/edu/duke/cs/osprey/spackstar/RunSpackStarDesign.java \
    src/test/java/edu/duke/cs/osprey/spackstar/RunWmbSequenceSubtreeValidation.java

JARGS="--add-opens java.base/java.util=ALL-UNNAMED \
       --add-opens java.base/java.lang=ALL-UNNAMED \
       --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
       -Xms4g -Xmx96g"

java $JARGS \
    -Dwmbis.validation.state="$STATE" \
    -Dwmbis.validation.numCpus="$SLURM_CPUS_PER_TASK" \
    -Dwmbis.validation.numGpus=0 \
    -Dwmbis.validation.iBounds="$IBOUNDS" \
    -Dwmbis.validation.samples="$SAMPLES" \
    -Dwmbis.validation.seed="$SEED" \
    -Dwmbis.validation.confidence="$CONFIDENCE" \
    -Dwmbis.validation.assignments="$ASSIGNMENTS" \
    -Dspackstar.wmbis.localCapMaxAssignments="$LOCAL_CAP_MAX_ASSIGNMENTS" \
    -Dspackstar.wmbis.probTailQuantile="$PROB_TAIL_QUANTILE" \
    -Dspackstar.run.mutableResidues="$MUTABLE_RESIDUES" \
    -Dspackstar.run.proteinWtFlexResidues="$PROTEIN_WT_FLEX" \
    -Dspackstar.run.ligandWtFlexResidues="$LIGAND_WT_FLEX" \
    -cp "$CP" \
    edu.duke.cs.osprey.spackstar.RunWmbSequenceSubtreeValidation

echo ""
echo "End: $(date)"
