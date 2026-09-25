#!/bin/bash
#SBATCH --job-name=build_markstar_million5
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
repo=/home/users/lz280/OSPREY3-fresh-packstar
root=/usr/xtmp/lz280/markstar_million5_20260918/build_$SLURM_JOB_ID
census=/usr/xtmp/lz280/packstar_million_census_20260918/build_12633933
mkdir -p "$root/source" "$root/classes" "$root/tmp"
export TMPDIR="$root/tmp"
cp "$repo/src/test/java/edu/duke/cs/osprey/packstar/MarkStarMillionControl.java" "$root/source/"
cp "$repo/slurm/scripts/run_markstar_million5_20260918.slurm" "$root/source/run.slurm"
cp "${BASH_SOURCE[0]}" "$root/source/build.slurm"
bash -n "$root/source/run.slurm"
/home/users/lz280/java/jdk-17.0.2+8/bin/javac -J-Xmx2g -cp "$(cat "$census/classpath.txt")" \
 -d "$root/classes" "$root/source/MarkStarMillionControl.java"
sha256sum "$root/source/"* > "$root/source.sha256"
touch "$root/READY"
sbatch --parsable --account=grisman --partition=grisman --array=0-4 \
 --export=ALL,CONTROL_BUILD="$root" "$root/source/run.slurm" > "$root/array_job.txt"
cat "$root/array_job.txt"
