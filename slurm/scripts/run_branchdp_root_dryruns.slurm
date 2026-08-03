#!/bin/bash
#SBATCH --job-name=branchdp_root_dryrun
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/branchdp_root_dryrun_%j.out
#SBATCH --error=slurm/logs/branchdp_root_dryrun_%j.err

set -euo pipefail

REPO=/home/users/lz280/IdeaProjects/OSPREY3
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
CP_FILE=/usr/xtmp/lz280/slurm_logs/.classpath_3k3q_gputest.txt
PDB_ROOT=/usr/xtmp/lz280/dance_bench/pdbs_prepped
CACHE_ROOT=/usr/xtmp/lz280/bench_comparison/results/emat_cache
RUN_ROOT=/tmp/branchdp-root-dryrun-${SLURM_JOB_ID}
ROOT_HOST_BUDGET_BYTES=${ROOT_HOST_BUDGET_BYTES:-751619276800}
DECOMP_STRATEGY=${DECOMP_STRATEGY:-hicks}
EXACT_IMPROVE_MAX_MILLIS=${EXACT_IMPROVE_MAX_MILLIS:-120000}
WEIGHTED_HICKS_RESTARTS=${WEIGHTED_HICKS_RESTARTS:-24}
WEIGHTED_HICKS_RANDOM_MOVES=${WEIGHTED_HICKS_RANDOM_MOVES:-32}
DP_TABLE_MODE=${DP_TABLE_MODE:-auto}
DP_MMAP_THRESHOLD_BYTES=${DP_MMAP_THRESHOLD_BYTES:-137438953472}
BENCH_METHOD=${BENCH_METHOD:-packstar}
PROFILE_DESIGN=${PROFILE_DESIGN:-3bua}
PROFILE_STATE=${PROFILE_STATE:-complex}
PROFILE_SEQ_INDEX=${PROFILE_SEQ_INDEX:-21}

cd "$REPO"
PROJECT_CACHE=/tmp/osprey3-project-cache-${SLURM_JOB_ID}
SKIP_COMPILE=${SKIP_COMPILE:-0}

echo "Branch-DP production-shape root dry-runs"
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURMD_NODENAME}  Start: $(date)"
echo "decompStrategy=$DECOMP_STRATEGY rootHostBudgetBytes=$ROOT_HOST_BUDGET_BYTES exactImproveMaxMillis=$EXACT_IMPROVE_MAX_MILLIS"
echo "weightedHicksRestarts=$WEIGHTED_HICKS_RESTARTS weightedHicksRandomMoves=$WEIGHTED_HICKS_RANDOM_MOVES"
echo "dpTableMode=$DP_TABLE_MODE dpMmapThresholdBytes=$DP_MMAP_THRESHOLD_BYTES"
echo "benchMethod=$BENCH_METHOD profileDesign=$PROFILE_DESIGN profileState=$PROFILE_STATE profileSeqIndex=$PROFILE_SEQ_INDEX"

if [ "$SKIP_COMPILE" = "1" ]; then
    echo "SKIP_COMPILE=1; using the already-validated workspace classes"
else
    export GRADLE_USER_HOME=/tmp/osprey3-gradle-${SLURM_JOB_ID}
    ./gradlew --project-cache-dir "$PROJECT_CACHE" \
        compileJava --no-daemon --max-workers=8
fi

echo "git_head=$(git rev-parse HEAD)"
sha256sum \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuFullDP.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.java \
    src/main/java/edu/duke/cs/osprey/branchdp/DPTable.java \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuFullDP.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPGpuOutOfCore.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/DPTable.class \
    build/classes/java/main/edu/duke/cs/osprey/branchdp/MappedDPTable.class \
    build/resources/main/gpuKernels/cuda/dp.bin

mkdir -p "$RUN_ROOT/emat_cache"
cp -a "$CACHE_ROOT/3bua_fennario" "$RUN_ROOT/emat_cache/"
cp -a "$CACHE_ROOT/3k3q_fennario" "$RUN_ROOT/emat_cache/"

CP=$(<"$CP_FILE")
JARGS=(
    --add-opens java.base/java.util=ALL-UNNAMED
    --add-opens java.base/java.lang=ALL-UNNAMED
    --add-opens java.base/java.lang.invoke=ALL-UNNAMED
    -Xmx56g
    -Xms4g
    -XX:-UseSuperWord
    -Dbranchdp.dp.dryRun=true
    -Dbranchdp.dp.cache=false
    -Dpackstar.dp.cache=false
    -Dbranchdp.dp.tableMode="$DP_TABLE_MODE"
    -Dbranchdp.dp.mmap.thresholdBytes="$DP_MMAP_THRESHOLD_BYTES"
    -Dbranchdp.rootSplit=gpubytes
    -Dbranchdp.rootSplit.gpuBudgetBytes=23085449216
    # Dry-run performs no DP allocation; use a production-scale logical budget
    # so candidate feasibility is measured independently of this job's 56 GiB heap.
    -Dbranchdp.rootSplit.hostBudgetBytes="$ROOT_HOST_BUDGET_BYTES"
    -Dbranchdp.decomp.strategy="$DECOMP_STRATEGY"
    -Dbranchdp.decomp.exactImprove.maxMillis="$EXACT_IMPROVE_MAX_MILLIS"
    -Dbranchdp.decomp.weightedHicks.restarts="$WEIGHTED_HICKS_RESTARTS"
    -Dbranchdp.decomp.weightedHicks.randomMoves="$WEIGHTED_HICKS_RANDOM_MOVES"
    -Dbranchdp.cutoff.strategy=RESIDUAL_BUDGET
    -Dbranchdp.cutoff.residualBudget=0.5
    -Dosprey.bench.method="$BENCH_METHOD"
    -Dosprey.bench.numCPUs=8
    -Dosprey.bench.outputDir="$RUN_ROOT"
)
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

if [ "$BENCH_METHOD" = "dp_profile" ]; then
    PROFILE_ARGS=(
        -Dosprey.dpProfile.state="$PROFILE_STATE"
        -Dosprey.dpProfile.seqIndex="$PROFILE_SEQ_INDEX"
        -Dosprey.dpProfile.maxMut=1
        -Dosprey.dpProfile.cacheFamily=packstar
    )
    if [ "$PROFILE_DESIGN" = "3bua" ]; then
        echo "BEGIN 3bua profile $(date)"
        /usr/bin/time -f 'ROOT_DRYRUN_WALL design=3bua-profile elapsed=%e maxRssKiB=%M' \
            "$JAVA" "${JARGS[@]}" "${PROFILE_ARGS[@]}" \
            -Dosprey.bench.designId=3bua_fennario \
            -Dosprey.bench.pdbPath="$PDB_ROOT/3bua/3bua.min.reduce.renum.pdb" \
            '-Dosprey.bench.mutable=C447;C450' \
            '-Dosprey.bench.flexible=D638;D641;D642;D645;D623;D619;D616;D620;C417;C414;C410;C421;C452;C443;C444;C446' \
            -cp "$CP" "$MAIN"
    elif [ "$PROFILE_DESIGN" = "3k3q" ]; then
        echo "BEGIN 3k3q profile $(date)"
        /usr/bin/time -f 'ROOT_DRYRUN_WALL design=3k3q-profile elapsed=%e maxRssKiB=%M' \
            "$JAVA" "${JARGS[@]}" "${PROFILE_ARGS[@]}" \
            -Dosprey.bench.designId=3k3q_fennario \
            -Dosprey.bench.pdbPath="$PDB_ROOT/3k3q/3k3q.min.reduce.renum.pdb" \
            '-Dosprey.bench.mutable=C444;C440' \
            '-Dosprey.bench.flexible=B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439' \
            -cp "$CP" "$MAIN"
    else
        echo "Unknown PROFILE_DESIGN=$PROFILE_DESIGN (expected 3bua or 3k3q)" >&2
        exit 2
    fi
    echo "Done: $(date)"
    exit 0
fi

echo "BEGIN 3bua $(date)"
/usr/bin/time -f 'ROOT_DRYRUN_WALL design=3bua elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" \
    -Dosprey.bench.designId=3bua_fennario \
    -Dosprey.bench.pdbPath="$PDB_ROOT/3bua/3bua.min.reduce.renum.pdb" \
    '-Dosprey.bench.mutable=C447;C450' \
    '-Dosprey.bench.flexible=D638;D641;D642;D645;D623;D619;D616;D620;C417;C414;C410;C421;C452;C443;C444;C446' \
    -cp "$CP" "$MAIN"
echo "END 3bua $(date)"

echo "BEGIN 3k3q $(date)"
/usr/bin/time -f 'ROOT_DRYRUN_WALL design=3k3q elapsed=%e maxRssKiB=%M' \
    "$JAVA" "${JARGS[@]}" \
    -Dosprey.bench.designId=3k3q_fennario \
    -Dosprey.bench.pdbPath="$PDB_ROOT/3k3q/3k3q.min.reduce.renum.pdb" \
    '-Dosprey.bench.mutable=C444;C440' \
    '-Dosprey.bench.flexible=B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439' \
    -cp "$CP" "$MAIN"
echo "END 3k3q $(date)"

echo "Done: $(date)"
