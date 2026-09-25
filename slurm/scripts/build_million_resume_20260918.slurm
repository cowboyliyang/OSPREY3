#!/bin/bash
#SBATCH --job-name=ps_ccd_resume_build
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
OLD=/usr/xtmp/lz280/packstar_million_census_20260918/build_12633933
ROOT=/usr/xtmp/lz280/packstar_million_census_20260918/resume_build_$SLURM_JOB_ID
REPO=/home/users/lz280/OSPREY3-fresh-packstar
mkdir -p "$ROOT/source" "$ROOT/classes" "$ROOT/tmp"
cp -a "$OLD/input" "$OLD/screen" "$ROOT/"
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/PackStarMillionCensus.java" "$ROOT/source/"
cp "$REPO/slurm/scripts/run_million_census_20260918.slurm" "$ROOT/source/"
cp "$OLD/input.sha256" "$ROOT/previous_inputs.sha256"
CP=$(cat "$OLD/classpath.txt")
JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin
"$JAVA/javac" -J-Xmx1g -cp "$CP" -d "$ROOT/classes" "$ROOT/source/PackStarMillionCensus.java"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/classpath.txt"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
ARGS=(--add-modules jdk.incubator.foreign --add-opens java.base/java.util=ALL-UNNAMED
 --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED
 -XX:-UseSuperWord -Xmx1g -Dcensus.cpus=2 "-Djava.io.tmpdir=$TMPDIR"
 -cp "$(cat "$ROOT/classpath.txt")" edu.duke.cs.osprey.packstar.PackStarMillionCensus)
"$JAVA/java" "${ARGS[@]}" shard "$ROOT" "$ROOT/check_original" "$ROOT/screen/4wyu.tsv" 0 8
"$JAVA/java" "${ARGS[@]}" shard "$ROOT" "$ROOT/check_resume" "$ROOT/screen/4wyu.tsv" 0 8 "$ROOT/check_original"
python3 - "$ROOT" <<'PY'
import csv, pathlib, struct, sys
r=pathlib.Path(sys.argv[1])
def load(p):
    data=(p/'energies.bin').read_bytes()
    assert len(data)==8*16
    return dict(struct.iter_unpack('>qd',data))
assert load(r/'check_original')==load(r/'check_resume')
row=next(csv.DictReader((r/'check_resume/resume.tsv').open(),delimiter='\t'))
assert int(row['restored'])==8
p=r/'check_partial';p.mkdir()
for name in ('layout.tsv','config.tsv'):(p/name).write_bytes((r/'check_original'/name).read_bytes())
(p/'energies.bin').write_bytes((r/'check_original/energies.bin').read_bytes()[:-5])
PY
"$JAVA/java" "${ARGS[@]}" shard "$ROOT" "$ROOT/check_tail" "$ROOT/screen/4wyu.tsv" 0 8 "$ROOT/check_partial"
python3 - "$ROOT" <<'PY'
import csv, pathlib, struct, sys
r=pathlib.Path(sys.argv[1])
def load(name):return dict(struct.iter_unpack('>qd',(r/name/'energies.bin').read_bytes()))
assert load('check_original')==load('check_tail')
row=next(csv.DictReader((r/'check_tail/resume.tsv').open(),delimiter='\t'))
assert int(row['restored'])==7 and int(row['trailing_bytes_recomputed'])==11
PY
sha256sum "$ROOT/source/"* > "$ROOT/source.sha256"
touch "$ROOT/READY"
echo "RESUME_BUILD=$ROOT"
