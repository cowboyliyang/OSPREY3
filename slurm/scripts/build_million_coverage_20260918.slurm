#!/bin/bash
#SBATCH --job-name=ps_mcov_build
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
CENSUS=/usr/xtmp/lz280/packstar_million_census_20260918/build_12633933
ROOT=/usr/xtmp/lz280/packstar_million_coverage_20260918/build_$SLURM_JOB_ID
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/tmp"
export TMPDIR="$ROOT/tmp"
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/PackStarMillionCoverage.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/"*million_coverage_20260918* "$ROOT/source/"
cp /usr/xtmp/lz280/packstar_small_z_20260918/build_12633897/profile.jvmargs "$ROOT/profile.jvmargs"
CP=$(cat "$CENSUS/classpath.txt")
/home/users/lz280/java/jdk-17.0.2+8/bin/javac -J-Xmx4g -cp "$CP" -d "$ROOT/classes" "$ROOT/source/PackStarMillionCoverage.java"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/classpath.txt"
sha256sum "$ROOT/source/PackStarMillionCoverage.java" "$ROOT/profile.jvmargs" "$CENSUS/input/"*.pdb "$CENSUS/screen/"*.tsv > "$ROOT/input.sha256"
touch "$ROOT/READY"
