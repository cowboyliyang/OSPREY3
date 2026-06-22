#!/bin/bash
#SBATCH --job-name=resBudget
#SBATCH --output=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_%j.out
#SBATCH --error=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -u

echo "=============================================="
echo "BranchMARK* residual-budget sweep (epsilon=0.68)"
echo "=============================================="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "Start:    $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.68}
BUDGETS=${BUDGETS:-"0.5 1.5 3.0 5.0"}

# ===== Clean Gradle locks/daemons =====
echo ""
echo "Cleaning Gradle locks and daemons..."
./gradlew --stop 2>/dev/null || true
pkill -u "$USER" -f "GradleDaemon" 2>/dev/null || true
rm -rf ~/.gradle/daemon/* 2>/dev/null || true
find .gradle -name "*.lock" -delete 2>/dev/null || true

# ===== Compile once =====
echo ""
echo "Compiling..."
./gradlew compileJava compileTestJava --no-daemon 2>&1 | tail -10
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "COMPILATION FAILED"
    exit 1
fi
echo "Compilation successful."

# ===== Sweep =====
SUMMARY=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_${SLURM_JOB_ID}_summary.tsv
printf "budget\tmarkstar_ms\tbranch_ms\tspeedup\tms_lb\tms_ub\tbr_lb\tbr_ub\n" > "$SUMMARY"

for BUDGET in $BUDGETS; do
    echo ""
    echo "=============================================="
    echo " residualBudget = $BUDGET   (numFlex=$NUM_FLEXIBLE, eps=$EPSILON)"
    echo "=============================================="

    LOG=/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/resBudget_${SLURM_JOB_ID}_b${BUDGET}.log

    TEST_JAVA_TOOL_OPTIONS=" \
      -Dbranchdp.test.numFlexible=$NUM_FLEXIBLE \
      -Dbranchdp.test.epsilon=$EPSILON \
      -Dbranchdp.cutoff.strategy=residual_budget \
      -Dbranchdp.cutoff.residualBudget=$BUDGET \
      -Dbranchdp.cutoff.keepConnected=true"

    JAVA_TOOL_OPTIONS="$TEST_JAVA_TOOL_OPTIONS" \
    ./gradlew test \
        --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testFlatSumVsOriginalMARKStar" \
        --no-daemon --rerun-tasks \
        -Dorg.gradle.jvmargs="-Xmx48g" \
        -Dbranchdp.test.numFlexible="$NUM_FLEXIBLE" \
        -Dbranchdp.test.epsilon="$EPSILON" \
        -Dbranchdp.cutoff.strategy=residual_budget \
        -Dbranchdp.cutoff.residualBudget="$BUDGET" \
        -Dbranchdp.cutoff.keepConnected=true \
        2>&1 | tee "$LOG"

    # Parse runtime/scores from the gradle stdout itself
    MS=$(grep -oE 'Original MARK\* +runtime: *[0-9]+ ms' "$LOG" | tail -1 | grep -oE '[0-9]+')
    BR=$(grep -oE 'BranchMARK\* FLAT_SUM +runtime: *[0-9]+ ms' "$LOG" | tail -1 | grep -oE '[0-9]+')
    SP=$(grep -oE 'Speedup \(MARK\*/FlatSum\): *[0-9.]+x' "$LOG" | tail -1 | grep -oE '[0-9.]+')
    MS_LB=$(grep -oE 'Original MARK\* +K\*: \[[^,]+' "$LOG" | tail -1 | grep -oE '[-0-9.]+$')
    MS_UB=$(grep -oE 'Original MARK\* +K\*: \[[^]]+\]' "$LOG" | tail -1 | grep -oE '[-0-9.]+\]' | tr -d ']')
    BR_LB=$(grep -oE 'BranchMARK\* FLAT_SUM +K\*: \[[^,]+' "$LOG" | tail -1 | grep -oE '[-0-9.]+$')
    BR_UB=$(grep -oE 'BranchMARK\* FLAT_SUM +K\*: \[[^]]+\]' "$LOG" | tail -1 | grep -oE '[-0-9.]+\]' | tr -d ']')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$BUDGET" "${MS:-NA}" "${BR:-NA}" "${SP:-NA}" \
        "${MS_LB:-NA}" "${MS_UB:-NA}" "${BR_LB:-NA}" "${BR_UB:-NA}" >> "$SUMMARY"
done

echo ""
echo "=============================================="
echo "SWEEP SUMMARY  (epsilon=$EPSILON, numFlex=$NUM_FLEXIBLE)"
echo "=============================================="
column -t "$SUMMARY"
echo ""
echo "End: $(date)"
echo "=============================================="
