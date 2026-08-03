#!/bin/bash
#SBATCH --job-name=wmbis_small
#SBATCH --output=/usr/xtmp/lz280/osprey3_slurm/logs/wmbis_small_%j.out
#SBATCH --error=/usr/xtmp/lz280/osprey3_slurm/logs/wmbis_small_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lz280@duke.edu

set -euo pipefail

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

echo "=============================================="
echo "WMB-IS small validation"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "Java: $(java -version 2>&1 | head -n 1)"
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

javac -cp "$CP" -d build/classes/java/main \
    src/main/java/edu/duke/cs/osprey/wmb/WeightedMiniBucket.java \
    src/main/java/edu/duke/cs/osprey/wmb/WmbModel.java

javac -cp "$CP" -d build/classes/java/test \
    src/test/java/edu/duke/cs/osprey/wmb/TestWeightedMiniBucket.java \
    src/test/java/edu/duke/cs/osprey/wmb/RunWmbSmallValidation.java

JARGS="--add-opens java.base/java.util=ALL-UNNAMED \
       --add-opens java.base/java.lang=ALL-UNNAMED \
       --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
       -Xms512m -Xmx4g"

java $JARGS -cp "$CP" edu.duke.cs.osprey.wmb.RunWmbSmallValidation

echo ""
echo "End: $(date)"
