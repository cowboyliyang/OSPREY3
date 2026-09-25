#!/bin/bash
#SBATCH --job-name=ps_z_build
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=fennario-03
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056
ROOT=/usr/xtmp/lz280/packstar_small_z_20260918/build_${SLURM_JOB_ID}
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/input" "$ROOT/tmp"
export TMPDIR="$ROOT/tmp"
# Reuse verified production bytecode; compile only the new experiment harness.
diff -qr "$REPO/src/main" "$BASE/source/src/main"
cp -a "$BASE/source/src" "$ROOT/source/"
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/PackStarSmallCoverage.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/"*small_z_20260918* "$ROOT/source/"
cp /home/users/lz280/BranchMARK\*_paper/recomb_draft/recomb_fs_v1.jvmargs "$ROOT/profile.jvmargs"
cp "$REPO/docs/PACKSTAR_SMALL_Z_COVERAGE_20260918.md" "$ROOT/protocol.md"
for system in 1gwc 2p4a 2rl0 4wyu; do
  cp "/usr/xtmp/lz280/dance_bench/pdbs_prepped/$system/$system.min.reduce.renum.pdb" "$ROOT/input/$system.pdb"
done
CP=$(cat "$BASE/test_classpath.txt")
/home/users/lz280/java/jdk-17.0.2+8/bin/javac -J-Xmx8g -cp "$CP" -d "$ROOT/classes" "$ROOT/source/PackStarSmallCoverage.java"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/classpath.txt"
sha256sum "$ROOT/source/PackStarSmallCoverage.java" "$ROOT/profile.jvmargs" "$ROOT/input/"*.pdb > "$ROOT/inputs.sha256"
cp "$BASE/source.sha256" "$ROOT/production_source.sha256"
git -C "$REPO" rev-parse HEAD > "$ROOT/git_head.txt"
touch "$ROOT/READY"
echo "BUILD_ROOT=$ROOT"
