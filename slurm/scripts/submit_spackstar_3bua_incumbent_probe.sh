#!/bin/bash
# Freeze, compile, checksum, and submit the 3BUA PACK* incumbent probe.

set -euo pipefail

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
SNAPSHOT_PARENT=/home/users/lz280/tmp/spackstar_source_snapshots
SLURM_SCRIPT="$REPO/slurm/scripts/run_spackstar_3bua_incumbent_probe.slurm"

mkdir -p "$SNAPSHOT_PARENT"
SNAPSHOT=$(mktemp -d "$SNAPSHOT_PARENT/3bua_incumbent.XXXXXXXX")
mkdir -p "$SNAPSHOT/classes"

cd "$REPO"
mapfile -t SOURCES < <(find src/main/java/edu/duke/cs/osprey/spackstar -name '*.java' -print | sort)
SOURCES+=(
    src/test/java/edu/duke/cs/osprey/spackstar/RunSpackStarDesign.java
    src/test/java/edu/duke/cs/osprey/spackstar/RunSpackStarEnvelopeProbe.java
    src/test/java/edu/duke/cs/osprey/spackstar/RunSpackStarIncumbentProbe.java
)

for source in "${SOURCES[@]}"; do
    mkdir -p "$SNAPSHOT/$(dirname "$source")"
    cp -- "$source" "$SNAPSHOT/$source"
done

cd "$SNAPSHOT"
find src -type f -name '*.java' -print0 \
    | sort -z \
    | xargs -0 sha256sum > source-manifest.sha256
SOURCE_DIGEST=$(sha256sum source-manifest.sha256 | awk '{print $1}')
RAW_CP=$(<"$CP_FILE")
IFS=':' read -r -a CP_ENTRIES <<< "$RAW_CP"
CP=
for entry in "${CP_ENTRIES[@]}"; do
    if [[ "$entry" != /* ]]; then
        entry="$REPO/$entry"
    fi
    if [[ -n "$CP" ]]; then
        CP+=:
    fi
    CP+="$entry"
done
printf '%s\n' "$CP" > classpath.txt

"$JAVA_HOME/bin/javac" -cp "$CP" -d classes "${SOURCES[@]}"

echo "$SOURCE_DIGEST" > source.digest
echo "Frozen SPACK* snapshot: $SNAPSHOT"
echo "SPACK* source digest:    $SOURCE_DIGEST"

sbatch \
    --export="ALL,SPACKSTAR_SNAPSHOT=$SNAPSHOT,SPACKSTAR_SOURCE_DIGEST=$SOURCE_DIGEST" \
    "$SLURM_SCRIPT"
