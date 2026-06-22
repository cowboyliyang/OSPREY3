#!/bin/bash
#SBATCH --job-name=resB
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_b%x_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_b%x_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --exclusive
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -u

BUDGET=${BUDGET:?must set BUDGET}
NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}

echo "=============================================="
echo "BranchMARK* residualBudget=$BUDGET (eps=$EPSILON, numFlex=$NUM_FLEXIBLE)"
echo "Job:  $SLURM_JOB_ID   Node: $SLURM_NODELIST   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

# Isolate BOTH ~/.gradle (GRADLE_USER_HOME) AND project-local .gradle (project-cache-dir)
# so the 6 jobs don't fight over fileHashes.lock
export GRADLE_USER_HOME="/tmp/lz280-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE="/tmp/lz280-projcache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE"

echo "Gradle user home:   $GRADLE_USER_HOME"
echo "Project cache dir:  $PROJECT_CACHE"
echo ""

# Compile was already done on the login node before submit, but redo per-job
# into the isolated cache so test execution doesn't replay it
echo "Compiling (isolated)..."
./gradlew --project-cache-dir "$PROJECT_CACHE" compileJava compileTestJava --no-daemon 2>&1 | tail -10
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "COMPILATION FAILED"
    rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
    exit 1
fi
echo "Compilation OK."
echo ""

TEST_JAVA_TOOL_OPTIONS=" \
  -Dbranchdp.test.numFlexible=$NUM_FLEXIBLE \
  -Dbranchdp.test.epsilon=$EPSILON \
  -Dbranchdp.cutoff.strategy=residual_budget \
  -Dbranchdp.cutoff.residualBudget=$BUDGET \
  -Dbranchdp.cutoff.keepConnected=true"

JAVA_TOOL_OPTIONS="$TEST_JAVA_TOOL_OPTIONS" \
./gradlew --project-cache-dir "$PROJECT_CACHE" test \
    --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testFlatSumVsOriginalMARKStar" \
    --no-daemon --rerun-tasks \
    -Dorg.gradle.jvmargs="-Xmx72g" \
    -Dbranchdp.test.numFlexible="$NUM_FLEXIBLE" \
    -Dbranchdp.test.epsilon="$EPSILON" \
    -Dbranchdp.cutoff.strategy=residual_budget \
    -Dbranchdp.cutoff.residualBudget="$BUDGET" \
    -Dbranchdp.cutoff.keepConnected=true \
    2>&1

echo ""
echo "End: $(date)"
rm -rf "$GRADLE_USER_HOME" "$PROJECT_CACHE"
