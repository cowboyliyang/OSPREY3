#!/bin/bash
#SBATCH --job-name=proposal_learning_test
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/proposal_learning_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/proposal_learning_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman
set -euo pipefail
cd /home/users/lz280/OSPREY3-fresh-packstar
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
run_root=/usr/xtmp/lz280/packstar_migration/proposal-learning-$SLURM_JOB_ID
mkdir -p "$run_root/tmp"
export TMPDIR="$run_root/tmp"
export OSPREY_GRADLE_XTMP_BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/debug-gradle-build
export GRADLE_USER_HOME=/usr/xtmp/lz280/packstar_migration/debug-gradle-home
./gradlew test --no-daemon --max-workers=1 \
    --project-cache-dir /usr/xtmp/lz280/packstar_migration/residue-verify-12505265/project-cache \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    --init-script slurm/scripts/verify_ccd.init.gradle \
    -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
    --tests edu.duke.cs.osprey.packstar.TestPackStarProposalLearning \
    --tests edu.duke.cs.osprey.packstar.TestPackStarTripleEtaCorrections \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC \
    --tests edu.duke.cs.osprey.packstar.TestPackStarConfig
cp -a "$OSPREY_GRADLE_XTMP_BUILD_ROOT"/*/test-results "$run_root/test-results"
