#!/bin/bash
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --nodelist=fennario-01
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --job-name=gnn_dbg2
#SBATCH --output=gnn_debug2_%j.out
#SBATCH --mail-user=lz280@duke.edu
#SBATCH --mail-type=END,FAIL

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo "---"

# Forward system properties to test JVM
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.benchmarkGNNStrategies" \
    -DtestMaxHeap=32g \
    --info \
    -Dorg.gradle.jvmargs="-Xmx1g" 2>&1 | head -5

# Use direct java command instead
CLASSPATH=$(./gradlew -q printTestClasspath 2>/dev/null || echo "")

# Build classpath manually
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do
    CP="$CP:$jar"
done

echo "Running test with strategy6..."
/home/users/lz280/java/jdk-17.0.2+8/bin/java \
    --add-opens java.base/java.util=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
    -Xmx32g -Xms2g \
    -Dosprey.gnn.benchStrategy=strategy6 \
    -Dosprey.gnn.numSeqs=1 \
    -cp "$CP" \
    org.junit.platform.console.ConsoleLauncher \
    --select-method "edu.duke.cs.osprey.markstar.TestBranchMARKStar#benchmarkGNNStrategies" \
    2>&1

echo "---"
echo "Done: $(date)"
