#!/bin/bash
# ===================================================================
# PACK* paper Fig 1 ("sample-size convergence") + Table 3 ("no-eta ablation").
#
# Fixed n=10 (branchdp.test.numFlexible), sweeps n_s = packstar.pac.samples over
# NS_LIST, and for each n_s runs BOTH strict PACK* (eta correction on) and the
# no-eta ablation (packstar.pac.etaEnabled=false, PACK*(no-eta) = q_m * mean(phi)
# baseline) so the resulting CSVs directly give Fig 1's two curves.
# Table 3 is just the n_s=1000 slice of the same data -- rerun with NS_LIST=1000
# alone if you only want that table point.
#
# PACK*-only (no MARK* leg -- MARK*'s epsilon doesn't depend on n_s/eta).
# One job per (n_s, eta) pair, each writing its own CSV to avoid append races:
#   results/sample_size_ns<ns>_eta<0|1>.csv
# Merge for plotting:
#   python3 -c "import pandas as pd,glob; \
#     pd.concat([pd.read_csv(f) for f in glob.glob('$OUTDIR/results/sample_size_ns*.csv')]) \
#       .to_csv('$OUTDIR/results/sample_size_all.csv', index=False)"
#
# Usage: bash bench_sample_size_convergence.sh ["100 200 500 1000 2000"] ["true false"]
#   env overrides: NUM_FLEXIBLE=10  EPSILON=0.683
#                  PACKSTAR_CONFIDENCE=0.05  PACKSTAR_RESIDUAL_BOUND=1.0  PACKSTAR_MAX_EST_SAMPLES=4000
#                  GRISMAN_CPUS=104  GRISMAN_MEM=400G  GRISMAN_GPUS=8
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=/usr/xtmp/lz280/slurm_logs
NS_LIST=${1:-"100 200 500 1000 2000"}
ETA_LIST=${2:-"true false"}

NUM_FLEXIBLE=${NUM_FLEXIBLE:-10}
EPSILON=${EPSILON:-0.683}
PACKSTAR_CONFIDENCE=${PACKSTAR_CONFIDENCE:-0.05}
PACKSTAR_RESIDUAL_BOUND=${PACKSTAR_RESIDUAL_BOUND:-1.0}
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
echo "$CP" > "$LOGDIR/.classpath_sample_size.txt"

# IMPORTANT: the estimation-round sample count actually used is N*-sizing's pick,
# capped by packstar.pac.maxEstSamples -- which defaults to a FIXED 4000 regardless
# of n_s, and N*-sizing stops as soon as it predicts hitting packstar.pac.targetEpsilon
# (default 0.683, the same loose target used everywhere else). At n_s<=2000 the
# estimator reaches that loose target with far fewer samples than n_s allows, so
# without the two overrides below, epsilon is flat vs. n_s regardless of the sweep
# (verified empirically). We force the FULL n_s-derived budget to be used by (a)
# capping maxEstSamples/unreachableCap at n_2 = n_s - floor(0.5 n_s) - floor(0.1 n_s)
# (paper's n_2 <= n_s - n_1 - n_0) and (b) setting targetEpsilon far tighter than
# reachable, so N*-sizing always falls back to the cap instead of stopping early.
PACKSTAR_TARGET_EPSILON_TIGHT=${PACKSTAR_TARGET_EPSILON_TIGHT:-0.001}

for NS in $NS_LIST; do
    N1=$(( NS * 5 / 10 )); N0=$(( NS / 10 )); ESTCAP=$(( NS - N1 - N0 ))
    for ETA in $ETA_LIST; do
        ETATAG=$([ "$ETA" = "true" ] && echo 1 || echo 0)
        OUTPUT_CSV="$OUTDIR/results/sample_size_ns${NS}_eta${ETATAG}.csv"
        JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
            --cpus-per-task="$GRISMAN_CPUS" --mem="$GRISMAN_MEM" --gres=gpu:a5000:"$GRISMAN_GPUS" \
            --time=14-00:00:00 --mail-user=lz280@duke.edu --mail-type=END,FAIL \
            --job-name=samplesize_ns${NS}_eta${ETATAG} \
            --output="$LOGDIR/samplesize_ns${NS}_eta${ETATAG}_%j.out" \
            --error="$LOGDIR/samplesize_ns${NS}_eta${ETATAG}_%j.err" \
            --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-$GRISMAN_CPUS}; $JAVA $JARGS \
                -Dbranchdp.test.numFlexible=$NUM_FLEXIBLE \
                -Dosprey.scalingpac.method=packstar \
                -Dosprey.scalingpac.epsilon=$EPSILON \
                -Dosprey.scalingpac.outputCsv=$OUTPUT_CSV \
                -Dosprey.scalingpac.numCPUs=\$RUN_CPUS \
                -Dosprey.branchdp.numCpus=\$RUN_CPUS \
                -Dpackstar.pac.samples=$NS \
                -Dpackstar.pac.confidence=$PACKSTAR_CONFIDENCE \
                -Dpackstar.pac.residualBound=$PACKSTAR_RESIDUAL_BOUND \
                -Dpackstar.pac.etaEnabled=$ETA \
                -Dpackstar.pac.maxEstSamples=$ESTCAP \
                -Dpackstar.pac.unreachableCap=$ESTCAP \
                -Dpackstar.pac.targetEpsilon=$PACKSTAR_TARGET_EPSILON_TIGHT \
                -cp \"\$(cat $LOGDIR/.classpath_sample_size.txt)\" $MAIN scaling_pac 2>&1" \
            | awk '{print $4}')
        echo "Submitted sample-size n_s=$NS eta=$ETA estCap=$ESTCAP: job $JID -> $OUTPUT_CSV"
    done
done
