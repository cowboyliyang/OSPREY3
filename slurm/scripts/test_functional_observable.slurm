#!/bin/bash
# Targeted Slurm validation for the PACK* functional observable path.
# All compilation and test execution stays off the login node.
#SBATCH --job-name=packstar_functional_observable
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/functional_observable_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/functional_observable_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=grisman
#SBATCH --account=grisman

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=/usr/xtmp/lz280/packstar_migration/functional-observable-gradle-build
GRADLE_ROOT=/usr/xtmp/lz280/packstar_migration/functional-observable-gradle-home

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_ROOT"
export GRADLE_USER_HOME="$GRADLE_ROOT"

cd "$REPO"
./gradlew compileJava compileTestJava \
    --no-daemon \
    --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g

./gradlew test \
    --no-daemon \
    --max-workers=1 \
    --init-script slurm/scripts/gradle_xtmp.init.gradle \
    -Dorg.gradle.jvmargs=-Xmx8g \
    --tests edu.duke.cs.osprey.packstar.TestPackStarFunctionalObservable
