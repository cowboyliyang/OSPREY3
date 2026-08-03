#!/bin/bash
# ===================================================================
# PACK* paper Table 2 ("scaling with n"): MARK* vs PACK* epsilon/CCD-count
# scaling on ConfSpaces2RL0.buildWildTypeConfSpace(n), n in {8,10,12,16,20}.
#
# Driver: PackStarScalingBench (src/test/java/.../markstar/bench/), invoked via
# TestBranchMARKStar main-dispatch mode "scaling_pac". One SLURM job per n value
# (each job runs BOTH the MARK* leg and the PACK* leg unless METHOD overrides),
# each job writes its OWN csv (results/scaling_n_<n>.csv) to avoid concurrent-
# append races between the 5 parallel jobs; merge afterwards for plotting:
#   python3 -c "import pandas as pd,glob; \
#     pd.concat([pd.read_csv(f) for f in glob.glob('$OUTDIR/results/scaling_n_*.csv')]) \
#       .to_csv('$OUTDIR/results/scaling_n_all.csv', index=False)"
#
# Usage: bash bench_scaling_pac.sh ["8 10 12 16 20"]
#   env overrides: EPSILON=0.683  METHOD=both|markstar|packstar
#                  PACKSTAR_SAMPLES=1000  PACKSTAR_CONFIDENCE=0.05  PACKSTAR_RESIDUAL_BOUND=1.0
#                  PACKSTAR_ETA_ENABLED=true  PACKSTAR_MAX_EST_SAMPLES=4000
#                  GRISMAN_CPUS=104  GRISMAN_MEM=400G  GRISMAN_GPUS=8
#                  JAVA_XMX=192g JAVA_XMS=192g
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=/usr/xtmp/lz280/slurm_logs
N_LIST=${1:-"8 10 12 16 20"}

EPSILON=${EPSILON:-0.683}
METHOD=${METHOD:-both}
PACKSTAR_SAMPLES=${PACKSTAR_SAMPLES:-1000}
PACKSTAR_CONFIDENCE=${PACKSTAR_CONFIDENCE:-0.05}
PACKSTAR_RESIDUAL_BOUND=${PACKSTAR_RESIDUAL_BOUND:-1.0}
PACKSTAR_ETA_ENABLED=${PACKSTAR_ETA_ENABLED:-true}
PACKSTAR_MAX_EST_SAMPLES=${PACKSTAR_MAX_EST_SAMPLES:-4000}
GRISMAN_CPUS=${GRISMAN_CPUS:-104}
GRISMAN_MEM=${GRISMAN_MEM:-400G}
GRISMAN_GPUS=${GRISMAN_GPUS:-8}
JAVA_XMX=${JAVA_XMX:-192g}
JAVA_XMS=${JAVA_XMS:-192g}

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
EXTRA_JVM_ARGS="${EXTRA_JVM_ARGS:-} -XX:-UseSuperWord -XX:+UseParallelGC -Dpackstar.dp.gpu=true -Dpackstar.dp.gpu.multiGpu=true -Dpackstar.dp.gpu.persistentContext=true"
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx$JAVA_XMX -Xms$JAVA_XMS $EXTRA_JVM_ARGS"
MAIN=edu.duke.cs.osprey.markstar.TestBranchMARKStar

mkdir -p "$LOGDIR" "$OUTDIR/results"

echo "Compiling test classes..."
./gradlew testClasses --no-daemon -Dorg.gradle.vfs.watch=false 2>&1 | tail -3
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do CP="$CP:$jar"; done
echo "$CP" > "$LOGDIR/.classpath_scaling_pac.txt"

# Same fix as bench_sample_size_convergence.sh: N*-sizing stops as soon as it
# predicts hitting the loose 0.683 target, so without pinning maxEstSamples/
# unreachableCap to n_s's own n_2 budget (and forcing the target unreachable),
# PACK*'s reported epsilon here would just reflect the 0.683 target rather than
# actually exercising the fixed n_s=$PACKSTAR_SAMPLES sample budget across n.
N1=$(( PACKSTAR_SAMPLES * 5 / 10 )); N0=$(( PACKSTAR_SAMPLES / 10 )); ESTCAP=$(( PACKSTAR_SAMPLES - N1 - N0 ))
PACKSTAR_TARGET_EPSILON_TIGHT=${PACKSTAR_TARGET_EPSILON_TIGHT:-0.001}

for N in $N_LIST; do
    OUTPUT_CSV="$OUTDIR/results/scaling_n_${N}.csv"
    JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
        --cpus-per-task="$GRISMAN_CPUS" --mem="$GRISMAN_MEM" --gres=gpu:a5000:"$GRISMAN_GPUS" \
        --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL \
        --job-name=scalingpac_n${N} \
        --output="$LOGDIR/scalingpac_n${N}_%j.out" --error="$LOGDIR/scalingpac_n${N}_%j.err" \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-$GRISMAN_CPUS}; $JAVA $JARGS \
            -Dbranchdp.test.numFlexible=$N \
            -Dosprey.scalingpac.method=$METHOD \
            -Dosprey.scalingpac.epsilon=$EPSILON \
            -Dosprey.scalingpac.outputCsv=$OUTPUT_CSV \
            -Dosprey.scalingpac.numCPUs=\$RUN_CPUS \
            -Dosprey.branchdp.numCpus=\$RUN_CPUS \
            -Dpackstar.pac.samples=$PACKSTAR_SAMPLES \
            -Dpackstar.pac.confidence=$PACKSTAR_CONFIDENCE \
            -Dpackstar.pac.residualBound=$PACKSTAR_RESIDUAL_BOUND \
            -Dpackstar.pac.etaEnabled=$PACKSTAR_ETA_ENABLED \
            -Dpackstar.pac.maxEstSamples=$ESTCAP \
            -Dpackstar.pac.unreachableCap=$ESTCAP \
            -Dpackstar.pac.targetEpsilon=$PACKSTAR_TARGET_EPSILON_TIGHT \
            -cp \"\$(cat $LOGDIR/.classpath_scaling_pac.txt)\" $MAIN scaling_pac 2>&1" \
        | awk '{print $4}')
    echo "Submitted scaling_pac n=$N (packstar estCap=$ESTCAP): job $JID -> $OUTPUT_CSV"
done
