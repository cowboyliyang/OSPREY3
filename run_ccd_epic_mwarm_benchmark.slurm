#!/bin/bash
#SBATCH --job-name=ccd_epic_mwarm
#SBATCH --output=ccd_epic_mwarm_%j.out
#SBATCH --error=ccd_epic_mwarm_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "=============================================="
echo "CCD vs EPIC vs Multi Warm Start Benchmark"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Clean Gradle locks and daemons
echo ""
echo "Cleaning Gradle locks and daemons..."
./gradlew --stop 2>/dev/null || true
pkill -u $USER -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

# Compile
echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1
if [ $? -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation successful."

# Run benchmark
echo ""
echo "=========================================="
echo "Running benchmarkCCDvsEPICvsMultiWarm"
echo "=========================================="
echo ""
./gradlew test \
    --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmarkCCDvsEPICvsMultiWarm" \
    --no-daemon \
    -Dorg.gradle.jvmargs="-Xmx28g" \
    --info 2>&1 | grep -v "^$" | grep -v "^w:" | grep -v "^>" | grep -v "^Note:" | grep -v "^warning:" | grep -v "Deprecated"

# Extract results from XML report
echo ""
echo "=========================================="
echo "Results from XML report:"
echo "=========================================="
XML=build/test-results/test/TEST-edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.xml
if [ -f "$XML" ]; then
    # Extract system-out content
    python3 -c "
import xml.etree.ElementTree as ET
tree = ET.parse('$XML')
root = tree.getroot()
for so in root.iter('system-out'):
    if so.text:
        print(so.text)
" 2>/dev/null || grep -oP '(?<=<!\[CDATA\[)[\s\S]*?(?=\]\]>)' "$XML"
else
    echo "No XML report found"
fi

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
