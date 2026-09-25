#!/bin/bash
#SBATCH --job-name=ps_moment_cpu
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${BUILD_ROOT:?}"
ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/moment_${SLURM_JOB_ID}
mkdir -p "$ROOT/classes" "$ROOT/tmp"
export TMPDIR="$ROOT/tmp"
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/MomentFitSpeedBench.java "$ROOT/"
CP=$(cat "$BUILD_ROOT/test_classpath.txt")
"$JAVA_ROOT/bin/javac" -cp "$CP" -d "$ROOT/classes" "$ROOT/MomentFitSpeedBench.java"
"$JAVA_ROOT/bin/java" -XX:ActiveProcessorCount=1 -Djava.io.tmpdir="$TMPDIR" \
    -cp "$ROOT/classes:$CP" edu.duke.cs.osprey.packstar.MomentFitSpeedBench > "$ROOT/result.tsv"
cat "$ROOT/result.tsv"
