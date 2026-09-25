#!/bin/bash
#SBATCH --job-name=ps_ablation_build
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --gres=gpu:titan_v:1
#SBATCH --exclude=jerry3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_forward_20260913/12588099
BUILD_ROOT=/usr/xtmp/lz280/packstar_proposal_ablation_20260917/build_${SLURM_JOB_ID:?}
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
mkdir -p "$BUILD_ROOT/source" "$BUILD_ROOT/classes" "$BUILD_ROOT/tmp" "$BUILD_ROOT/runtime"
export TMPDIR="$BUILD_ROOT/tmp"
export PYTHONPYCACHEPREFIX="$BUILD_ROOT/pycache"

# About 60 MiB / a few thousand source files. Reuse only verified unchanged
# Kotlin and native/CUDA dependencies; every Java source is compiled anew.
rsync -a --exclude=/buildSrc/build/ --exclude=/buildSrc/.gradle/ \
    "$REPO/src" "$REPO/lib" "$REPO/buildSrc" "$REPO/gradle" "$REPO/slurm" \
    "$REPO/gradlew" "$REPO/build.gradle.kts" "$REPO/settings.gradle.kts" "$REPO/docs" \
    "$BUILD_ROOT/source/"
git -C "$REPO" rev-parse HEAD > "$BUILD_ROOT/git_head.txt"
git -C "$REPO" status --short > "$BUILD_ROOT/git_status.txt"
git -C "$REPO" diff --binary > "$BUILD_ROOT/working_tree.patch"
cd "$BUILD_ROOT/source"
diff -qr src/main/kotlin "$BASE/source/src/main/kotlin"
diff -qr src/test/kotlin "$BASE/source/src/test/kotlin"
diff -qr src/main/cc "$BASE/source/src/main/cc"
diff -qr src/main/cu "$BASE/source/src/main/cu"
while IFS= read -r kernel; do
    cmp "$kernel" "$BASE/source/$kernel"
done < <(rg --files src/main/resources/gpuKernels -g '*.cu')
for kind in main test; do
    mkdir -p "$BUILD_ROOT/runtime/$kind-resources" "$BUILD_ROOT/runtime/$kind-kotlin"
    rsync -a "$BASE/build/osprey3-root/resources/$kind/" "$BUILD_ROOT/runtime/$kind-resources/"
    rsync -a --exclude='*.bin' "src/$kind/resources/" "$BUILD_ROOT/runtime/$kind-resources/"
    rsync -a "$BASE/build/osprey3-root/classes/kotlin/$kind/" "$BUILD_ROOT/runtime/$kind-kotlin/"
done
CP="$BUILD_ROOT/runtime/main-kotlin:$BUILD_ROOT/runtime/test-kotlin:$BUILD_ROOT/runtime/main-resources:$BUILD_ROOT/runtime/test-resources"
IFS=: read -ra BASE_ENTRIES < "$BASE/test_classpath.txt"
for entry in "${BASE_ENTRIES[@]}"; do
    if [[ "$entry" == *.jar ]]; then
        test -s "$entry"
        CP="$CP:$entry"
    fi
done
rg --files src/main/java src/test/java -g '*.java' > "$BUILD_ROOT/java-sources.txt"
"$JAVA_ROOT/bin/javac" -J-Xmx20g --add-modules jdk.incubator.foreign \
    -cp "$CP" -d "$BUILD_ROOT/classes" @"$BUILD_ROOT/java-sources.txt"
printf '%s:%s\n' "$BUILD_ROOT/classes" "$CP" > "$BUILD_ROOT/test_classpath.txt"
"$JAVA_ROOT/bin/javac" -cp "$CP" -d "$BUILD_ROOT/classes" slurm/scripts/ForwardTripleTestRunner.java
"$JAVA_ROOT/bin/java" --add-modules jdk.incubator.foreign -Xmx16g \
    -Djava.io.tmpdir="$TMPDIR" -cp "$(cat "$BUILD_ROOT/test_classpath.txt")" ForwardTripleTestRunner
/bin/bash -n slurm/scripts/run_proposal_ablation_titan_20260917.slurm
/bin/bash -n slurm/scripts/gate_proposal_ablation_titan_20260917.slurm
/bin/bash -n slurm/scripts/summarize_proposal_ablation_titan_20260917.slurm
python3 -m py_compile slurm/scripts/summarize_proposal_ablation_titan_20260917.py slurm/scripts/check_proposal_ablation_20260917.py
rg --files src slurm/scripts -0 | sort -z | xargs -0 sha256sum > "$BUILD_ROOT/source.sha256"
{
    printf 'buildJob\t%s\n' "$SLURM_JOB_ID"
    printf 'baseDependencies\t%s\n' "$BASE"
    printf 'precision\tFP64\n'
    printf 'javaSources\tall-current-main-and-test\n'
    printf 'finishedAt\t%s\n' "$(date -Is)"
} > "$BUILD_ROOT/build_manifest.tsv"
touch "$BUILD_ROOT/READY"
echo "BUILD_READY=$BUILD_ROOT"
