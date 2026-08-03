#!/bin/bash
# Slurm-only compile, unit-test, and shell-contract validation for the
# sequence-bundle PACK* preflight implementation.
#
# This job intentionally does not launch a production preflight.  After it
# succeeds, submit run_packstar_preflight_multinode.slurm with CLASS_ROOT set
# to the xtmp build output printed at the end of this job.
#
#SBATCH --job-name=packstar_mn_check
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
INIT="$REPO/slurm/scripts/gradle_xtmp.init.gradle"
GRADLE_USER_HOME=${GRADLE_USER_HOME:-/usr/xtmp/lz280/gradle-home}
PACKSTAR_BUILD_ROOT=${PACKSTAR_BUILD_ROOT:-/usr/xtmp/lz280/osprey3-gradle-build/packstar_${SLURM_JOB_ID}}
OSPREY_GRADLE_XTMP_BUILD_ROOT=$PACKSTAR_BUILD_ROOT
export GRADLE_USER_HOME OSPREY_GRADLE_XTMP_BUILD_ROOT

cd "$REPO"
mkdir -p /usr/xtmp/lz280/slurm_logs "$GRADLE_USER_HOME"
PACKSTAR_TEST_ROOT=/usr/xtmp/lz280/packstar_validation/$SLURM_JOB_ID
mkdir -p "$PACKSTAR_TEST_ROOT/tmp"
export TMPDIR="$PACKSTAR_TEST_ROOT/tmp"

# Syntax checks are executed inside this Slurm job, in accordance with the
# workspace policy that project logic must not run on the login node.
bash -n slurm/scripts/packstar_preflight_sequence_rank.sh
bash -n slurm/scripts/aggregate_packstar_preflight_sequence.sh
bash -n slurm/scripts/merge_packstar_multinode_results.sh
bash -n slurm/scripts/run_packstar_preflight_multinode.slurm
bash -n slurm/scripts/run_packstar_multigpu_multinode.slurm
bash -n slurm/scripts/run_packstar_single_gpu_node.slurm
bash -n slurm/scripts/run_packstar_multinode_production.slurm
bash -n slurm/scripts/recover_packstar_multinode_results.slurm
bash -n slurm/scripts/compare_packstar_merged.slurm

python3 -m unittest discover -s slurm/wrappers/tests \
    -p 'test_packstar_*.py'

./gradlew --no-daemon --init-script "$INIT" \
    test \
    --tests edu.duke.cs.osprey.branchdp.TestBranchDpAdmission \
    --tests edu.duke.cs.osprey.packstar.TestPackStarAdmissionDecision \
    --tests edu.duke.cs.osprey.packstar.TestPackStarRandomStreams \
    --tests edu.duke.cs.osprey.kstar.TestKStarSequenceSharding

echo "compiledClassRoot=$PACKSTAR_BUILD_ROOT/osprey-root"
echo "packstar multinode validation passed job=$SLURM_JOB_ID time=$(date --iso-8601=seconds)"
