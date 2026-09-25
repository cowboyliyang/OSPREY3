#!/bin/bash
#SBATCH --job-name=ps_million_build
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056
ROOT=/usr/xtmp/lz280/packstar_million_census_20260918/build_$SLURM_JOB_ID
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/input" "$ROOT/tmp" "$ROOT/screen"
export TMPDIR="$ROOT/tmp"
diff -qr "$REPO/src/main" "$BASE/source/src/main"
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/PackStarMillionCensus.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/"*million_census_20260918* "$ROOT/source/"
for system in 1gwc 2p4a 2rl0 4wyu; do
 cp "/usr/xtmp/lz280/dance_bench/pdbs_prepped/$system/$system.min.reduce.renum.pdb" "$ROOT/input/$system.pdb"
done
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
CP=$(cat "$BASE/test_classpath.txt")
"$JAVA_ROOT/bin/javac" -J-Xmx4g -cp "$CP" -d "$ROOT/classes" "$ROOT/source/PackStarMillionCensus.java"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/classpath.txt"
sha256sum "$ROOT/source/PackStarMillionCensus.java" "$ROOT/input/"*.pdb > "$ROOT/input.sha256"
cp "$BASE/source.sha256" "$ROOT/production_source.sha256"
"$JAVA_ROOT/bin/java" -Xmx6g -Djava.io.tmpdir="$TMPDIR" -cp "$(cat "$ROOT/classpath.txt")" \
 edu.duke.cs.osprey.packstar.PackStarMillionCensus screen "$ROOT" "$ROOT/screen"
touch "$ROOT/READY"
echo "BUILD_ROOT=$ROOT"
