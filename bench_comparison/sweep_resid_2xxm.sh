#!/bin/bash
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3
OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
CP=$(cat $LOGDIR/.classpath_bench.txt)
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
PDB=/usr/xtmp/lz280/dance_bench/pdbs_prepped/2xxm/2xxm.min.reduce.renum.pdb
MUT='B107'
FLEX='B115;B117;B109;B173;B105;B167;B165;A30;A4;A6;A5'
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench
PAC_SAMPLES=1000
Q=0.85
DPGPU="-Dbranchmarkstar.dp.gpu=true -Dbranchmarkstar.dp.gpu.multiGpu=true -Dbranchmarkstar.dp.gpu.persistentContext=true"
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx192g -Xms4g -XX:-UseSuperWord $DPGPU"

for R in 0.3 0.4 0.5 0.7 1.0 1.5; do
  RTAG=${R/./}
  RES=$OUTDIR/results/residsweep/r$RTAG
  mkdir -p "$RES"
  JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
        --gres=gpu:a5000:8 --cpus-per-task=104 --mem=128G \
        --time=06:00:00 --mail-type=FAIL --mail-user=lz280@duke.edu \
        --job-name=resid${RTAG}_2xxm \
        --output=$LOGDIR/pac_2xxm_resid${RTAG}_%j.out \
        --error=$LOGDIR/pac_2xxm_resid${RTAG}_%j.err \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-104}; $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$PDB \
            -Dosprey.bench.mutable='$MUT' \
            -Dosprey.bench.flexible='$FLEX' \
            -Dosprey.bench.method=pac \
            -Dosprey.bench.designId=2xxm \
            -Dosprey.bench.outputDir=$RES \
            -Dosprey.bench.numCPUs=\$RUN_CPUS \
            -Dbranchmarkstar.usePAC=true \
            -Dbranchmarkstar.pac.residualBound=$R \
            -Dbranchmarkstar.rootSplit=memory \
            -Dbranchmarkstar.dp.cache=false \
            -Dbranchmarkstar.pac.samples=$PAC_SAMPLES \
            -Dbranchmarkstar.pac.confidence=0.05 \
            -Dbranchmarkstar.pac.clip=true \
            -Dbranchmarkstar.pac.clipQuantile=$Q \
            -cp \"$CP\" $MAIN 2>&1" | awk '{print $4}')
  echo "submitted residualBound=$R -> job $JID  (results: $RES)"
done
