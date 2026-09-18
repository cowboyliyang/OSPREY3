#!/bin/bash
# Supply account, partition, GPU allocation and log paths to sbatch.
#SBATCH --job-name=packstar-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${REPO:?Absolute source checkout}"
: "${BUILD_ROOT:?New build directory on scratch}"
: "${JAVA_HOME:?JDK 17 installation}"
test ! -e "$BUILD_ROOT"
mkdir -p "$BUILD_ROOT/source" "$BUILD_ROOT/tmp"
export TMPDIR="$BUILD_ROOT/tmp"
export PYTHONPYCACHEPREFIX="$BUILD_ROOT/pycache"
export CUDA_CACHE_PATH="$BUILD_ROOT/cuda-cache"
export GRADLE_USER_HOME="${GRADLE_USER_HOME:-$BUILD_ROOT/gradle-home}"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT/build"
export PACKSTAR_CLASSPATH_FILE="$BUILD_ROOT/test_classpath.txt"
export OSPREY_CUDA_ARCHS="${OSPREY_CUDA_ARCHS:-90}"
export PATH="$JAVA_HOME/bin:$PATH"
# Source snapshot: roughly 100 MiB; dependency caches can require several GiB.
git -C "$REPO" rev-parse HEAD > "$BUILD_ROOT/git_head.txt"
git -C "$REPO" diff --binary HEAD > "$BUILD_ROOT/source.patch"
git -C "$REPO" ls-files -z | tar -C "$REPO" --null -T - -cf - | tar -C "$BUILD_ROOT/source" -xf -
cd "$BUILD_ROOT/source"
nvidia-smi > "$BUILD_ROOT/gpus.txt"
for file in slurm/h200/*.slurm; do bash -n "$file"; done
python3 -m py_compile slurm/h200/run.py
./gradlew compileCuda --no-daemon --max-workers=1 \
    -I slurm/scripts/gradle_xtmp.init.gradle -Dorg.gradle.jvmargs=-Xmx8g
./gradlew test --no-daemon --max-workers=1 \
    -I slurm/scripts/gradle_xtmp.init.gradle -I slurm/h200/classpath.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g -DtestMaxHeap=16g \
    --tests 'edu.duke.cs.osprey.packstar.TestPackStar*' \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies \
    --tests edu.duke.cs.osprey.branchdp.TestHigherOrderFullDP \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_nonLeafRegularMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_childSlicedMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_hybridMatchesJava \
    --tests edu.duke.cs.osprey.branchdp.TestGpuFullDP.gpuExactTriples_outOfCoreMatchesJava \
    writeExperimentClasspath
python3 - <<'PY'
import os
from pathlib import Path
import xml.etree.ElementTree as ET
reports = list(Path(os.environ['OSPREY_GRADLE_XTMP_BUILD_ROOT']).rglob('TEST-*.xml'))
assert reports, 'Missing JUnit reports'
count = 0
gpu = 0
for report in reports:
    suite = ET.parse(report).getroot()
    for case in suite.findall('testcase'):
        assert case.find('failure') is None and case.find('error') is None, ET.tostring(case)
        if 'TestGpuFullDP' in case.get('classname', ''):
            assert case.find('skipped') is None, 'GPU verification was skipped'
            gpu += 1
        count += 1
assert gpu == 4, f'Expected four GPU route checks; found {gpu}'
print(f'VERIFIED {count} tests, including {gpu} GPU route checks')
PY
find src slurm/h200 -type f -print0 | sort -z | xargs -0 sha256sum > "$BUILD_ROOT/source.sha256"
touch "$BUILD_ROOT/READY"
echo "BUILD_READY=$BUILD_ROOT"
