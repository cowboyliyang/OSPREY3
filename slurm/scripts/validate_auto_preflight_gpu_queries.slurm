#!/bin/bash
# One changed Java class plus regression runner; outputs under 1 MiB plus test logs.
#SBATCH --job-name=ps-auto-query-fix
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_auto_preflight/check_12705173
ROOT=/usr/xtmp/lz280/packstar_auto_preflight/query_fix_${SLURM_JOB_ID}
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/tmp"
export TMPDIR="$ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
CP=$(<"$BASE/test_classpath.txt")
cp "$REPO/src/main/java/edu/duke/cs/osprey/branchdp/BranchDpBackend.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/AutoPreflightRegressionGate.java" "$ROOT/source/"
javac -cp "$CP" -d "$ROOT/classes" "$ROOT/source/BranchDpBackend.java" "$ROOT/source/AutoPreflightRegressionGate.java"
java --add-modules jdk.incubator.foreign -Xmx8g -Djava.io.tmpdir="$TMPDIR" \
    -cp "$ROOT/classes:$CP" AutoPreflightRegressionGate
printf '%s\n' "$ROOT/classes:$CP" > "$ROOT/test_classpath.txt"
ln -sfn "$ROOT" /usr/xtmp/lz280/packstar_auto_preflight/current
echo "AUTO_PREFLIGHT_CLASSPATH=$ROOT/test_classpath.txt"
