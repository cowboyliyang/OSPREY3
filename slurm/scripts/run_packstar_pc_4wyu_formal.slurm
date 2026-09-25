#!/bin/bash
# Formal 4WYU A+D PACK* P_C design case, bounded protocol v3.
# All compilation, PDB splitting, minimization, enumeration, and analysis run
# inside this Slurm job. Home contains only source and configuration inputs.
#SBATCH --job-name=packstar_pc_4wyu
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=grisman
#SBATCH --account=grisman

set -euo pipefail
test -n "${SLURM_JOB_ID:-}"

REPO=/home/users/lz280/OSPREY3-fresh-packstar
SOURCE_PDB=/usr/xtmp/lz280/dance_bench/pdbs_prepped/4wyu/4wyu.min.reduce.renum.pdb
EXPECTED_SHA256=db330b3e3968865b8225c3f7735df44d39828c3bd032bc2f2ce3d3def7c77d53
RUN_ROOT=${RUN_ROOT:-/usr/xtmp/lz280/packstar_functional_validation/formal_4wyu_${SLURM_JOB_ID}}
if [ -n "${VALIDATED_RUN:-}" ]; then
    test -f "$VALIDATED_RUN/VALIDATED"
    exec bash "$VALIDATED_RUN/source_snapshot/slurm/scripts/launch_packstar_pc_4wyu.sh" "$VALIDATED_RUN" "$RUN_ROOT"
fi
BUILD_ROOT=$RUN_ROOT/build
GRADLE_CACHE_SOURCE=/usr/xtmp/lz280/packstar_rescue_20260905/gradle_cache
GRADLE_ROOT=${GRADLE_ROOT:-$RUN_ROOT/gradle_home}
JAVA_HOME=${JAVA_HOME:-/home/users/lz280/java/jdk-17.0.2+8}

mkdir -p "$RUN_ROOT"/input "$RUN_ROOT"/output \
    "$RUN_ROOT"/source_snapshot "$RUN_ROOT"/frequency_severity "$RUN_ROOT"/tmp
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME
export PATH="$JAVA_HOME/bin:$PATH"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT/gradle_build"
export GRADLE_USER_HOME="$GRADLE_ROOT"
mkdir -p "$GRADLE_ROOT"
if [ ! -d "$GRADLE_ROOT/caches/modules-2" ]; then
    mkdir -p "$GRADLE_ROOT/caches"
    cp -a "$GRADLE_CACHE_SOURCE/wrapper" "$GRADLE_ROOT/"
    cp -a "$GRADLE_CACHE_SOURCE/caches/modules-2" "$GRADLE_ROOT/caches/"
fi

cd "$REPO"
printf 'job_id\t%s\n' "$SLURM_JOB_ID" > "$RUN_ROOT/job_manifest.tsv"
printf 'node\t%s\n' "${SLURMD_NODENAME:-unknown}" >> "$RUN_ROOT/job_manifest.tsv"
printf 'started_at\t%s\n' "$(date --iso-8601=seconds)" >> "$RUN_ROOT/job_manifest.tsv"
printf 'source_pdb\t%s\n' "$SOURCE_PDB" >> "$RUN_ROOT/job_manifest.tsv"
printf 'expected_source_sha256\t%s\n' "$EXPECTED_SHA256" >> "$RUN_ROOT/job_manifest.tsv"

test -s "$SOURCE_PDB"
SOURCE_SHA256=$(sha256sum "$SOURCE_PDB" | awk '{print $1}')
printf 'actual_source_sha256\t%s\n' "$SOURCE_SHA256" >> "$RUN_ROOT/job_manifest.tsv"
test "$SOURCE_SHA256" = "$EXPECTED_SHA256"
cp docs/packstar_pc_4wyu_formal_v3.json "$RUN_ROOT/protocol.json"

# Preserve the exact dirty source used by the job and create the two explicit
# state inputs from the frozen prepared PDB. Chain B/C are rejected below.
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ --exclude=/.gradle/ --exclude='*.class' \
    src lib buildSrc gradle slurm gradlew build.gradle.kts settings.gradle.kts \
    "$RUN_ROOT/source_snapshot/"
git diff --binary > "$RUN_ROOT/tracked_changes.patch"
git status --short > "$RUN_ROOT/source_status.txt"
git rev-parse HEAD > "$RUN_ROOT/git_head.txt"
mkdir -p "$RUN_ROOT/source_snapshot/examples/python.KStar" "$RUN_ROOT/source_snapshot/docs"
cp examples/python.KStar/2RL0.min.reduce.pdb "$RUN_ROOT/source_snapshot/examples/python.KStar/"
cp docs/packstar_pc_4wyu_formal_v3.json "$RUN_ROOT/source_snapshot/docs/"
cd "$RUN_ROOT/source_snapshot"
find src slurm docs -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_ROOT/source_sha256.txt"

awk 'BEGIN { found=0 } /^(ATOM  |HETATM)/ && substr($0,22,1)=="A" { print; found=1 } END { if (!found) exit 1; print "END" }' \
    "$SOURCE_PDB" > "$RUN_ROOT/input/protein_A.pdb"
awk 'BEGIN { found=0 } /^(ATOM  |HETATM)/ && substr($0,22,1)=="D" { print; found=1 } END { if (!found) exit 1; print "END" }' \
    "$SOURCE_PDB" > "$RUN_ROOT/input/peptide_D.pdb"
test -s "$RUN_ROOT/input/protein_A.pdb"
test -s "$RUN_ROOT/input/peptide_D.pdb"
sha256sum "$RUN_ROOT/input/protein_A.pdb" "$RUN_ROOT/input/peptide_D.pdb" \
    > "$RUN_ROOT/input/input_sha256.txt"
awk '/^(ATOM  |HETATM)/ { chains[substr($0,22,1)]=1 } END { for (c in chains) print c }' \
    "$SOURCE_PDB" | sort > "$RUN_ROOT/input/source_chains.txt"
test "$(wc -l < "$RUN_ROOT/input/source_chains.txt")" -eq 4
test "$(cat "$RUN_ROOT/input/source_chains.txt")" = $'A\nB\nC\nD'

# Compile and obtain a runtime class path entirely in xtmp. The classpath
# task is registered by the repository's Slurm init script.
export PACKSTAR_FORMAL_CLASSPATH_FILE="$RUN_ROOT/classpath.txt"
./gradlew test --offline --no-daemon --max-workers=1 \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalGeometry \
    --tests edu.duke.cs.osprey.packstar.TestPackStarResult \
    writeFormalCaseClasspath \
    --project-cache-dir "$RUN_ROOT/project-cache" \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/packstar_formal_classpath.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -DtestMaxHeap=8g \
    -Dorg.gradle.jvmargs="-Xmx16g -Djava.io.tmpdir=$TMPDIR"
test -s "$RUN_ROOT/classpath.txt"
# writeFormalCaseClasspath already verifies compileTestJava's actual output.
# Do not assume the root project's name is the same in a source snapshot.

exec bash slurm/scripts/launch_packstar_pc_4wyu.sh "$RUN_ROOT" "$RUN_ROOT"
