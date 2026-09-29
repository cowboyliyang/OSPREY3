#!/bin/bash
#SBATCH --job-name=ps-push-check
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_push_20260929/check_${SLURM_JOB_ID}
# Source snapshot approximately 100 MiB / several thousand files; artifacts on xtmp.
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$RUN_ROOT/build"
export PACKSTAR_CLASSPATH_FILE="$RUN_ROOT/test_classpath.txt"
cd "$REPO"
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts "$RUN_ROOT/source/"
git diff --binary HEAD > "$RUN_ROOT/source.patch"
git ls-files --others --exclude-standard > "$RUN_ROOT/untracked_sources.txt"
cd "$RUN_ROOT/source"
python3 - "$RUN_ROOT" <<'PY'
import ast, pathlib, subprocess, sys
root=pathlib.Path(sys.argv[1])
paths=list((root/'source/slurm').rglob('*.slurm'))
for p in paths:subprocess.run(['bash','-n',str(p)],check=True)
for p in (root/'source/slurm').rglob('*.py'):ast.parse(p.read_text(),filename=str(p))
print('PASS: Slurm shell and Python source syntax')
PY
./gradlew test --offline --no-daemon --max-workers=1 \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    -I slurm/scripts/gradle_xtmp.init.gradle -I slurm/scripts/verify_ccd.init.gradle \
    -I slurm/h200/classpath.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" -DtestMaxHeap=8g \
    --tests edu.duke.cs.osprey.kstar.TestKStarScore \
    --tests edu.duke.cs.osprey.kstar.TestKStarScoreLog10 \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC \
    --tests edu.duke.cs.osprey.packstar.TestPackStarAutoPreflight \
    --tests edu.duke.cs.osprey.packstar.TestPackStarAdmissionDecision \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDpAdmission \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleDecompositionCosts \
    writeExperimentClasspath
echo "VERIFIED_PUSH_SOURCE=$RUN_ROOT"
