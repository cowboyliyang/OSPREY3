#!/bin/bash
#SBATCH --job-name=ccd3cal_build
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588
ORIGINAL=/usr/xtmp/lz280/packstar_pair38_20260921/launch_12677939
ROOT=/usr/xtmp/lz280/packstar_3cal_trp_protein_census_20260924/build_$SLURM_JOB_ID
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/input" "$ROOT/tmp" "$ROOT/pilot_sparse" "$ROOT/pilot_full"
export TMPDIR="$ROOT/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/PackStar3calProteinCensus.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/"*3cal_trp_census_20260924* "$ROOT/source/"
cp "$REPO/slurm/scripts/merge_million_census_20260918.py" "$ROOT/source/"
cp "$ORIGINAL/inputs/3cal/3cal.min.reduce.renum.pdb" "$ROOT/input/3cal.pdb"
cp "$ORIGINAL/runs/3cal_full_pair-only_s42_J12678012_T16/manifest.json" "$ROOT/input/original_manifest.json"
cp "$ORIGINAL/runs/3cal_full_pair-only_s42_J12678012_T16/adaptive_frequency_severity/state-e7ef13207bdd8a21/frequency_severity_final_samples.tsv" "$ROOT/input/original_samples.tsv"
cp "$BASE/source/src/test/java/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.java" "$ROOT/source/GenericPDBBench.original.java"
python3 - "$ROOT" <<'PY'
import hashlib,sys
from pathlib import Path
root=Path(sys.argv[1])
assert hashlib.sha256((root/'input/3cal.pdb').read_bytes()).hexdigest()=='ae9d613aad4f87b80cd00c5aa3dfcf518dbd91c3aa007f6f28ae191f264f37f1'
src=(root/'source/GenericPDBBench.original.java').read_text()
src=src.replace('public class GenericPDBBench {','public class Census3calBench {')
src=src.replace('public static void main(String[] args)', 'public static TestKStar.ConfSpaces buildSpaces(String[] args)')
start=src.index('        int wmbGpus = Integer.getInteger("osprey.wmb.numGpus", 0);')
end=src.index('    private static void exportGNNData(')
src=src[:start]+'        return confSpaces;\n    }\n\n'+src[end:]
(root/'source/Census3calBench.java').write_text(src)
PY
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
CP=$(cat "$BASE/test_classpath.txt")
"$JAVA_ROOT/bin/javac" -J-Xmx4g -cp "$CP" -d "$ROOT/classes" \
 "$ROOT/source/PackStar3calProteinCensus.java" "$ROOT/source/Census3calBench.java"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/classpath.txt"
cp "$BASE/source.sha256" "$ROOT/production_source.sha256"
ARGS=(--add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED
 -XX:-UseSuperWord -Xmx6g "-Djava.io.tmpdir=$TMPDIR" -cp "$(cat "$ROOT/classpath.txt")")
"$JAVA_ROOT/bin/java" "${ARGS[@]}" -Dosprey.bench.method=ccd_census \
 "-Dosprey.bench.pdbPath=$ROOT/input/3cal.pdb" -Dosprey.bench.mutable=C147 \
 '-Dosprey.bench.flexible=D210;D212;C126;C128;C113;C137;C111;C139;C120;C149' \
 -Dosprey.bench.designId=3cal "-Dosprey.bench.outputDir=$ROOT" -Dosprey.bench.numCPUs=4 \
 "-Dcensus.originalSamples=$ROOT/input/original_samples.tsv" \
 edu.duke.cs.osprey.packstar.PackStar3calProteinCensus prepare "$ROOT" "$ROOT" > "$ROOT/prepare.log" 2>&1
for TARGET in sparse full; do
 "$JAVA_ROOT/bin/java" "${ARGS[@]}" -Dcensus.cpus=4 "-Dcensus.target=$TARGET" \
  edu.duke.cs.osprey.packstar.PackStar3calProteinCensus pilot "$ROOT" "$ROOT/pilot_$TARGET" 0 256 > "$ROOT/pilot_$TARGET/run.log" 2>&1
done
sha256sum "$ROOT/source/"* "$ROOT/input/"* "$ROOT/eref.ser" > "$ROOT/inputs.sha256"
touch "$ROOT/READY"
cat "$ROOT/original_energy_check.txt" "$ROOT/pilot_sparse/result.tsv" "$ROOT/pilot_full/result.tsv"
echo "BUILD_ROOT=$ROOT"
