#!/bin/bash
#SBATCH --job-name=verify_residue_ccd
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/verify_residue_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/verify_residue_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
set -euo pipefail
cd /home/users/lz280/OSPREY3-fresh-packstar
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$JAVA_HOME/bin:$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
export NVCC="$cuda_env/bin/nvcc"
run_root=/usr/xtmp/lz280/packstar_migration/residue-verify-$SLURM_JOB_ID
mkdir -p "$run_root/tmp"
export TMPDIR="$run_root/tmp"
export OSPREY_GRADLE_XTMP_BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/debug-gradle-build
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/debug-gradle-home
gradle_args=(--no-daemon --max-workers=1
    --project-cache-dir "$run_root/project-cache"
    --init-script slurm/scripts/gradle_xtmp.init.gradle
    "-Dorg.gradle.jvmargs=-Xmx8g -Djava.io.tmpdir=$TMPDIR")
# Run after the debug compile job. Keep the normal resource-probing path on.
for count in 1 2 8; do
    ./gradlew "${gradle_args[@]}" -Dosprey.cuda.debugBatch=true \
        -Dosprey.cuda.test.batchCount="$count" test \
        --tests edu.duke.cs.osprey.minimization.TestMinimization.testResidueCudaCCDBatchMatchesScalar \
        --tests edu.duke.cs.osprey.minimization.TestMinimization.testResidueCudaCCDPartialBatchMatchesScalar
    cp -a "$OSPREY_GRADLE_XTMP_BUILD_ROOT"/*/test-results "$run_root/batch-$count-results"
done
./gradlew "${gradle_args[@]}" test \
    --tests edu.duke.cs.osprey.minimization.TestMinimization.auditResidueCudaCCDAgainstCpu \
    --tests edu.duke.cs.osprey.packstar.TestPackStarConformationDedup
cp -a "$OSPREY_GRADLE_XTMP_BUILD_ROOT"/*/test-results "$run_root/cpu-audit-results"
