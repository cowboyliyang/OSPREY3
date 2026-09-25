#!/bin/bash
#SBATCH --job-name=final_evalue_test
#SBATCH --account=grisman
#SBATCH --partition=compsci,grisman
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588
RUN_ROOT=/usr/xtmp/lz280/packstar_final_evalue_validation_20260924_$SLURM_JOB_ID
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
# Six source snapshots plus classes, logs and replay tables: below 20 MiB.
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/classes" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
cp "$REPO/src/main/java/edu/duke/cs/osprey/packstar/PackStarEstimator.java" "$RUN_ROOT/source/"
cp "$REPO/src/main/java/edu/duke/cs/osprey/packstar/PackStarFrequencySeverityPAC.java" "$RUN_ROOT/source/"
cp "$REPO/src/main/java/edu/duke/cs/osprey/energy/ResidueInteractions.java" "$RUN_ROOT/source/"
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/TestPackStarFrequencySeverityPAC.java" "$RUN_ROOT/source/"
cp "$REPO/slurm/scripts/FinalSeverityTestRunner.java" "$RUN_ROOT/source/"
cp "$REPO/slurm/scripts/audit_final_evalue_20260924.py" "$RUN_ROOT/source/"
CP=$(cat "$BASE/test_classpath.txt")
"$JAVA_ROOT/bin/javac" -J-Xmx4g -cp "$CP" -sourcepath '' -d "$RUN_ROOT/classes" \
    "$RUN_ROOT/source/ResidueInteractions.java" \
    "$RUN_ROOT/source/PackStarFrequencySeverityPAC.java" \
    "$RUN_ROOT/source/PackStarEstimator.java" \
    "$RUN_ROOT/source/TestPackStarFrequencySeverityPAC.java" \
    "$RUN_ROOT/source/FinalSeverityTestRunner.java"
python3 "$RUN_ROOT/source/audit_final_evalue_20260924.py" > "$RUN_ROOT/audit.log"
AUDIT_ROOT=/usr/xtmp/lz280/packstar_final_evalue_audit_20260924_$SLURM_JOB_ID
"$JAVA_ROOT/bin/java" --add-modules jdk.incubator.foreign \
    --add-opens java.base/java.util=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
    -Xmx4g "-Djava.io.tmpdir=$TMPDIR" -cp "$RUN_ROOT/classes:$CP" \
    FinalSeverityTestRunner "$AUDIT_ROOT/java_replay.tsv" > "$RUN_ROOT/tests.log" 2>&1
sha256sum "$RUN_ROOT/source/"* > "$RUN_ROOT/source.sha256"
cat "$RUN_ROOT/tests.log" "$AUDIT_ROOT/report.md"
echo "VALIDATION_ROOT=$RUN_ROOT"
echo "AUDIT_ROOT=$AUDIT_ROOT"
