#!/bin/bash
#SBATCH --job-name=ps_million_z
#SBATCH --account=grisman
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
CENSUS=${CENSUS_ROOT:-/usr/xtmp/lz280/packstar_million_census_20260918/build_12633933}
SYSTEM=${SYSTEM:-4wyu}
SEED=$((10000+${SLURM_ARRAY_TASK_ID:?}))
OUT="${RESULT_ROOT:?}/seed_$SEED"
mkdir -p "$OUT/tmp" "$OUT/cuda_cache"
test ! -e "$OUT/STARTED"
touch "$OUT/STARTED"
export TMPDIR="$OUT/tmp" CUDA_CACHE_PATH="$OUT/cuda_cache" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
nvidia-smi --query-gpu=index,name,uuid,memory.total,driver_version --format=csv > "$OUT/gpu.csv"
scontrol show job "$SLURM_JOB_ID" > "$OUT/slurm.txt"
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
ARGS=(--add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED
 -Xmx12g -XX:-UseSuperWord "@$BUILD_ROOT/profile.jvmargs" "-Djava.io.tmpdir=$TMPDIR"
 -Dbranchdp.dp.cache=false -Dpackstar.dp.cache=false
 -Dbranchdp.dp.gpu=true -Dbranchdp.dp.gpu.failIfNoGpuPath=true
 -Dbranchdp.dp.gpu.multiGpu=false -Dbranchdp.dp.gpu.maxGpus=1
 -Dbranchdp.dp.parallel.threads=4 -Dpackstar.dp.parallel.threads=4
 -Dbranchdp.rootSplit=gpubytes -Dbranchdp.rootSplit.hostBudgetBytes=8GiB -Dbranchdp.rootSplit.gpuBudgetBytes=6GiB
 -Dbranchdp.decomp.strategy=weighted_hicks -Dbranchdp.decomp.weightedHicks.restarts=24
 -Dbranchdp.decomp.weightedHicks.randomMoves=32
 -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET -Dbranchdp.cutoff.residualBudget=1
 -Dpackstar.pac.frequencySeverity.tripleEtaSelectionStrategy=budget-forward
 -Dpackstar.pac.frequencySeverity.tripleEtaMaxHostBytes=8GiB
 -Dpackstar.pac.sampling.threads=4 -Dpackstar.pac.sampling.gpu=false
 -Dpackstar.pac.targetEpsilon=0.683 "-Dpackstar.pac.randomSeed=$SEED"
 "-Dpackstar.pac.frequencySeverity.outputDir=$OUT/audit"
 -cp "$(cat "$BUILD_ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStarMillionCoverage "$CENSUS" "$OUT" "$SYSTEM")
printf '%q ' "$JAVA" "${ARGS[@]}" > "$OUT/command.sh"
set +e
/usr/bin/time -o "$OUT/wall.time" -f 'elapsed=%e maxRssKiB=%M' "$JAVA" "${ARGS[@]}" > "$OUT/run.log" 2>&1
RC=$?
set -e
printf '%s\n' "$RC" > "$OUT/exit_code"
exit "$RC"
