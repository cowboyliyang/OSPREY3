#!/bin/bash
#SBATCH --job-name=ps_z_coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
test -f "${BUILD_ROOT:?}/READY"
MODE=${MODE:-run}
TASK=${SLURM_ARRAY_TASK_ID:-0}
if [[ "$MODE" == reference ]]; then CASE=$TASK; SEED=0; else CASE=$((TASK / 100)); SEED=$((10000 + TASK % 100)); fi
SYSTEMS=(1gwc 2p4a 2rl0 4wyu)
SYSTEM=${SYSTEMS[$CASE]}
ROOT=${RESULT_ROOT:?}
if [[ "$MODE" == reference ]]; then OUT="$ROOT/reference/$SYSTEM"; else OUT="$ROOT/runs/$SYSTEM/seed_$SEED"; fi
mkdir -p "$OUT"
test ! -e "$OUT/STARTED"
touch "$OUT/STARTED"
mkdir -p "$OUT/tmp" "$OUT/cuda_cache"
export TMPDIR="$OUT/tmp" CUDA_CACHE_PATH="$OUT/cuda_cache" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
nvidia-smi --query-gpu=index,name,uuid,memory.total,driver_version --format=csv > "$OUT/gpu.csv"
scontrol show job "$SLURM_JOB_ID" > "$OUT/slurm.txt"
ARGS=(--add-modules jdk.incubator.foreign
 --add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED
 --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx12g -XX:-UseSuperWord
 "@$BUILD_ROOT/profile.jvmargs" "-Djava.io.tmpdir=$TMPDIR"
 "-Dcoverage.pdb=$BUILD_ROOT/input/$SYSTEM.pdb" -Dcoverage.cpus=4
 -Dbranchdp.dp.cache=false -Dpackstar.dp.cache=false
 -Dbranchdp.dp.gpu=true -Dbranchdp.dp.gpu.failIfNoGpuPath=true
 -Dbranchdp.dp.gpu.multiGpu=false -Dbranchdp.dp.gpu.maxGpus=1
 -Dbranchdp.dp.parallel.threads=4 -Dpackstar.dp.parallel.threads=4
 -Dbranchdp.rootSplit=gpubytes -Dbranchdp.rootSplit.hostBudgetBytes=8GiB
 -Dbranchdp.rootSplit.gpuBudgetBytes=6GiB
 -Dbranchdp.decomp.strategy=weighted_hicks -Dbranchdp.decomp.weightedHicks.restarts=24
 -Dbranchdp.decomp.weightedHicks.randomMoves=32
 -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET -Dbranchdp.cutoff.residualBudget=1
 -Dpackstar.pac.frequencySeverity.tripleEtaSelectionStrategy=budget-forward
 -Dpackstar.pac.frequencySeverity.tripleEtaMaxHostBytes=8GiB
 -Dpackstar.pac.sampling.threads=4 -Dpackstar.pac.sampling.gpu=false
 -Dpackstar.pac.targetEpsilon=0.683 "-Dpackstar.pac.randomSeed=$SEED"
 "-Dpackstar.pac.frequencySeverity.outputDir=$OUT/audit"
 -cp "$(cat "$BUILD_ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStarSmallCoverage
 "$MODE" "$CASE" "$OUT" "$ROOT/reference/$SYSTEM")
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
printf '%q ' "$JAVA" "${ARGS[@]}" > "$OUT/command.sh"
set +e
/usr/bin/time -o "$OUT/wall.time" -f 'elapsed=%e maxRssKiB=%M' "$JAVA" "${ARGS[@]}" > "$OUT/run.log" 2>&1
RC=$?
set -e
printf '%s\n' "$RC" > "$OUT/exit_code"
exit "$RC"
