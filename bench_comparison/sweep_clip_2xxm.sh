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

PACKSTAR_SAMPLES=${PACKSTAR_SAMPLES:-${PAC_SAMPLES:-1000}}
RESID=1.5
DPGPU="-Dpackstar.dp.gpu=true -Dpackstar.dp.gpu.multiGpu=true -Dpackstar.dp.gpu.persistentContext=true"
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx192g -Xms4g -XX:-UseSuperWord $DPGPU"

for Q in 0.80 0.85 0.90 0.93 0.95 0.97; do
  QTAG=${Q/./}
  RES=$OUTDIR/results/clipsweep/q$QTAG
  mkdir -p "$RES"
  JID=$(sbatch --partition=grisman --account=grisman --constraint=a5000 \
        --gres=gpu:a5000:8 --cpus-per-task=104 --mem=128G \
        --time=06:00:00 --mail-type=FAIL --mail-user=lz280@duke.edu \
        --job-name=clip${QTAG}_2xxm \
        --output=$LOGDIR/packstar_2xxm_clipq${QTAG}_%j.out \
        --error=$LOGDIR/packstar_2xxm_clipq${QTAG}_%j.err \
        --wrap "RUN_CPUS=\${SLURM_CPUS_ON_NODE:-104}; $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$PDB \
            -Dosprey.bench.mutable='$MUT' \
            -Dosprey.bench.flexible='$FLEX' \
            -Dosprey.bench.method=packstar \
            -Dosprey.bench.designId=2xxm \
            -Dosprey.bench.outputDir=$RES \
            -Dosprey.bench.numCPUs=\$RUN_CPUS \
            -Dpackstar.pac.residualBound=$RESID \
            -Dpackstar.rootSplit=memory \
            -Dpackstar.dp.cache=false \
            -Dpackstar.pac.samples=$PACKSTAR_SAMPLES \
            -Dpackstar.pac.confidence=0.05 \
            -Dpackstar.pac.clip=true \
            -Dpackstar.pac.clipQuantile=$Q \
            -cp \"$CP\" $MAIN 2>&1" | awk '{print $4}')
  echo "submitted clipQ=$Q -> job $JID  (results: $RES)"
done
