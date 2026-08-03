#!/bin/bash
# ===================================================================
# PACK* paper "PAC coverage verification": run PACK* 100x with different
# random seeds at n=10, n_s=1000, delta=0.05; report the fraction of runs
# whose interval [q-,q+] contains the ground-truth q (expected >=95%).
# Ground truth = MARK* run to a very low epsilon (or exhaustive enumeration
# on the sparse graph, not automated here) on the SAME conf space/seed-0 DP.
#
# Two steps, run in order:
#   bash bench_pac_coverage.sh ground_truth   # 1 job, MARK* at tight epsilon
#   bash bench_pac_coverage.sh seeds          # 100-task array job, PACK* per seed
#   bash bench_pac_coverage.sh all            # both (default)
# Then aggregate with: python3 analyze_pac_coverage.py
#
# Seed array jobs are sized small (1 GPU, ~1/8 node each) so up to
# GRISMAN_GPUS_TOTAL tasks run concurrently across the fennario a5000 fleet,
# instead of 100 serialized full-node jobs.
#
# env overrides: NUM_FLEXIBLE=10  PACKSTAR_SAMPLES=1000  PACKSTAR_CONFIDENCE=0.05
#                GROUND_TRUTH_EPSILON=0.02  PACKSTAR_RESIDUAL_BOUND=1.0
#                NUM_SEEDS=100  GRISMAN_GPUS_TOTAL=48
#                SEED_CPUS=12  SEED_MEM=45G
#                GRISMAN_CPUS=104  GRISMAN_MEM=400G  GRISMAN_GPUS=8   (ground_truth job only)
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=/usr/xtmp/lz280/slurm_logs
MODE=${1:-all}

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
PACKSTAR_SAMPLES=${PACKSTAR_SAMPLES:-1000}
PACKSTAR_CONFIDENCE=${PACKSTAR_CONFIDENCE:-0.05}
PACKSTAR_RESIDUAL_BOUND=${PACKSTAR_RESIDUAL_BOUND:-1.0}
GROUND_TRUTH_EPSILON=${GROUND_TRUTH_EPSILON:-0.02}
NUM_SEEDS=${NUM_SEEDS:-100}
GRISMAN_GPUS_TOTAL=${GRISMAN_GPUS_TOTAL:-48}
SEED_CPUS=${SEED_CPUS:-12}
SEED_MEM=${SEED_MEM:-45G}
GRISMAN_CPUS=${GRISMAN_CPUS:-104}
GRISMAN_MEM=${GRISMAN_MEM:-400G}
GRISMAN_GPUS=${GRISMAN_GPUS:-8}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-192g}
SEED_JAVA_XMX=${SEED_JAVA_XMX:-24g}
SEED_JAVA_XMS=${SEED_JAVA_XMS:-24g}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
MAIN=edu.duke.cs.osprey.markstar.TestBranchMARKStar

mkdir -p "$LOGDIR" "$OUTDIR/results"

echo "Compiling test classes..."
./gradlew testClasses --no-daemon -Dorg.gradle.vfs.watch=false 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_pac_coverage.txt"

submit_ground_truth() {
    local EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord -XX:+UseParallelGC"
    local JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx$JAVA_XMX -Xms$JAVA_XMS $EXTRA_JVM_ARGS"
    local OUTPUT_CSV="$OUTDIR/results/pac_coverage_ground_truth_n${NUM_FLEXIBLE}.csv"
    JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
        --cpus-per-task="$GRISMAN_CPUS" --mem="$GRISMAN_MEM" --gres=gpu:a5000:"$GRISMAN_GPUS" \
        --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL \
        --job-name=covGroundTruth \
        --output="$LOGDIR/covGroundTruth_%j.out" --error="$LOGDIR/covGroundTruth_%j.err" \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-$GRISMAN_CPUS}; $JAVA $JARGS \
            -Dbranchdp.test.numFlexible=$NUM_FLEXIBLE \
            -Dosprey.scalingpac.method=markstar \
            -Dosprey.scalingpac.epsilon=$GROUND_TRUTH_EPSILON \
            -Dosprey.scalingpac.outputCsv=$OUTPUT_CSV \
            -Dosprey.branchdp.numCpus=\$RUN_CPUS \
            -cp \"\$(cat $LOGDIR/.classpath_pac_coverage.txt)\" $MAIN scaling_pac 2>&1" \
        | awk '{print $4}')
    echo "Submitted coverage ground-truth (MARK* eps=$GROUND_TRUTH_EPSILON): job $JID -> $OUTPUT_CSV"
}

submit_seeds() {
    local EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord -XX:+UseParallelGC -Dpackstar.dp.gpu=true -Dpackstar.dp.gpu.multiGpu=false -Dpackstar.dp.gpu.persistentContext=true"
    local JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx$SEED_JAVA_XMX -Xms$SEED_JAVA_XMS $EXTRA_JVM_ARGS"
    # Same fix as bench_sample_size_convergence.sh: force the N*-sizing estimation
    # round to actually use the full n_s-derived budget (n_2 = n_s - n_1 - n_0)
    # instead of stopping early once it predicts hitting the loose 0.683 target.
    local N1=$(( PACKSTAR_SAMPLES * 5 / 10 )) N0=$(( PACKSTAR_SAMPLES / 10 ))
    local ESTCAP=$(( PACKSTAR_SAMPLES - N1 - N0 ))
    local TARGET_EPS_TIGHT=${PACKSTAR_TARGET_EPSILON_TIGHT:-0.001}
    JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
        --cpus-per-task="$SEED_CPUS" --mem="$SEED_MEM" --gres=gpu:a5000:1 \
        --array=1-${NUM_SEEDS}%${GRISMAN_GPUS_TOTAL} \
        --time=8:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL \
        --job-name=covSeed \
        --output="$LOGDIR/covSeed_%a_%A.out" --error="$LOGDIR/covSeed_%a_%A.err" \
        --wrap "SEED=\$SLURM_ARRAY_TASK_ID; OUTPUT_CSV=$OUTDIR/results/pac_coverage_seed_n${NUM_FLEXIBLE}_\${SEED}.csv; \
            $JAVA $JARGS \
            -Dbranchdp.test.numFlexible=$NUM_FLEXIBLE \
            -Dosprey.scalingpac.method=packstar \
            -Dosprey.scalingpac.epsilon=0.683 \
            -Dosprey.scalingpac.outputCsv=\$OUTPUT_CSV \
            -Dosprey.scalingpac.numCPUs=$SEED_CPUS \
            -Dpackstar.pac.samples=$PACKSTAR_SAMPLES \
            -Dpackstar.pac.confidence=$PACKSTAR_CONFIDENCE \
            -Dpackstar.pac.residualBound=$PACKSTAR_RESIDUAL_BOUND \
            -Dpackstar.pac.etaEnabled=true \
            -Dpackstar.pac.maxEstSamples=$ESTCAP \
            -Dpackstar.pac.unreachableCap=$ESTCAP \
            -Dpackstar.pac.targetEpsilon=$TARGET_EPS_TIGHT \
            -Dpackstar.pac.randomSeed=\$SEED \
            -cp \"\$(cat $LOGDIR/.classpath_pac_coverage.txt)\" $MAIN scaling_pac 2>&1" \
        | awk '{print $4}')
    echo "Submitted coverage seed array (1..$NUM_SEEDS, throttle $GRISMAN_GPUS_TOTAL, estCap=$ESTCAP): job $JID -> $OUTDIR/results/pac_coverage_seed_n${NUM_FLEXIBLE}_<seed>.csv"
}

case "$MODE" in
    ground_truth) submit_ground_truth ;;
    seeds)        submit_seeds ;;
    all)          submit_ground_truth; submit_seeds ;;
    *) echo "Usage: $0 [ground_truth|seeds|all]"; exit 1 ;;
esac
