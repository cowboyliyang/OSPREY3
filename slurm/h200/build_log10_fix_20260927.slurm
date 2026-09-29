#!/bin/bash
#SBATCH --job-name=ps-log10-fix-check
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_fix_check_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_fix_check_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
REPO=/home/users/lz280/OSPREY3-fresh-packstar
BASE=/usr/xtmp/lz280/packstar_pro6000_20260925/duke_12701382
ORIGINAL=/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588
export ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/build_${SLURM_JOB_ID}
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
mkdir -p "$ROOT"
mkdir "$ROOT/tmp" "$ROOT/classes" "$ROOT/test-classes" "$ROOT/source" "$ROOT/original_manifests"
test -f "$BASE/READY"
cp -a "$BASE/source/slurm" "$ROOT/source/"
mkdir -p "$ROOT/source/src/main/java/edu/duke/cs/osprey/kstar" "$ROOT/source/src/test/java/edu/duke/cs/osprey/kstar"
cp "$REPO/src/main/java/edu/duke/cs/osprey/kstar/KStarScore.java" "$ROOT/source/src/main/java/edu/duke/cs/osprey/kstar/"
cp "$REPO/src/test/java/edu/duke/cs/osprey/kstar/"{TestKStarScore.java,TestKStarScoreLog10.java} "$ROOT/source/src/test/java/edu/duke/cs/osprey/kstar/"
cp "$REPO/slurm/h200/"{Log10OutputRegression.java,select_log10_reruns_20260927.py,audit_all_log10_20260927.py,select_all_log10_reruns_20260927.py,node_log10_reruns_20260927.slurm} "$ROOT/source/slurm/h200/"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
cp "$BASE/git_head.txt" "$ROOT/git_head.txt"
cp "$BASE/source.sha256" "$ROOT/base-source.sha256"
cp "$ORIGINAL/source.sha256" "$ROOT/base-runtime-source.sha256"
cp "$BASE/test_classpath.txt" "$ROOT/base_classpath.txt"
# Use exactly the per-case script already copied by the running jobs.
cmp /usr/xtmp/lz280/packstar_pro6000_20260925/node_12701384/design.slurm /usr/xtmp/lz280/packstar_pro6000_20260925/node_12701385/design.slurm
cp /usr/xtmp/lz280/packstar_pro6000_20260925/node_12701385/design.slurm "$ROOT/design.slurm"
python3 - "$ORIGINAL" <<'PY'
import difflib, os, sys
from pathlib import Path
root = Path(os.environ['ROOT'])
relative = 'src/main/java/edu/duke/cs/osprey/kstar/KStarScore.java'
before = (Path(sys.argv[1]) / 'source' / relative).read_text()
after = (root / 'source' / relative).read_text()
start = '\tpublic static Double scoreToLog10(BigDecimal score) {'
end = '\tpublic static String scoreToLog10String(BigDecimal score) {'
def excluding_output_method(source):
    prefix, rest = source.split(start)
    _, suffix = rest.split(end)
    return prefix, suffix
assert excluding_output_method(before) == excluding_output_method(after)
(root / 'output-only.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True), fromfile='frozen/KStarScore.java', tofile='patched/KStarScore.java')))
print('PASS: only scoreToLog10 differs from the frozen production source')
PY
python3 "$ROOT/source/slurm/h200/select_log10_reruns_20260927.py" "$ROOT" | tee "$ROOT/selection.txt"
python3 - <<'PY'
import csv, hashlib, json, os, shutil
from pathlib import Path
root = Path(os.environ['ROOT'])
config = root / 'source/slurm/h200'
runner_sha = hashlib.sha256((config / 'run.py').read_bytes()).hexdigest()
manifest_sha = hashlib.sha256((config / 'frontier_active.tsv').read_bytes()).hexdigest()
expected = dict(cohort='frontier', mode='full', arm='pair-only', seed=42, rb=1, gpus=4, heap_gib=850, host_gib=800, gpu_gib=85)
production = dict(line.split('=', 1) for line in (config / 'production.properties').read_text().splitlines() if line and not line.startswith('#'))
production.update({'branchdp.cutoff.residualBudget': '1', 'packstar.pac.frequencySeverity.tripleEta': 'false'})
for row in csv.DictReader((root / 'cases.tsv').open(), delimiter='\t'):
    path = Path(row['original_path']) / 'manifest.json'
    m = json.loads(path.read_text())
    assert all(m['options'][k] == v for k, v in expected.items()), row
    assert m['cpus'] == 128 and m['precision'] == 'FP64' and m['ccd'] == 'CPU'
    assert m['runner_sha256'] == runner_sha and m['frontier_manifest_sha256'] == manifest_sha
    assert all(m['properties'][k] == v for k, v in production.items()), row
    shutil.copy2(path, root / 'original_manifests' / f"{row['index']}.json")
print('PASS: original runner, design/sequence manifest, all production properties and resource options match')
PY
CP=$(cat "$BASE/test_classpath.txt")
"$JAVA_HOME/bin/javac" -implicit:none -sourcepath "$ROOT/source/src/main/java" -cp "$CP" -d "$ROOT/classes" \
    "$ROOT/source/src/main/java/edu/duke/cs/osprey/kstar/KStarScore.java"
"$JAVA_HOME/bin/javac" -implicit:none -cp "$ROOT/classes:$CP" -d "$ROOT/test-classes" \
    "$ROOT/source/src/test/java/edu/duke/cs/osprey/kstar/"{TestKStarScore.java,TestKStarScoreLog10.java} \
    "$ROOT/source/slurm/h200/Log10OutputRegression.java"
timeout --kill-after=5s 120s "$JAVA_HOME/bin/java" --add-modules jdk.incubator.foreign -Xmx4g \
    -Djava.io.tmpdir="$TMPDIR" -cp "$ROOT/classes:$ROOT/test-classes:$CP" \
    Log10OutputRegression "$ROOT/archived_z.tsv" | tee "$ROOT/regression.txt"
test "$(find "$ROOT/classes" -name '*.class' | wc -l)" -eq 1
test -f "$ROOT/classes/edu/duke/cs/osprey/kstar/KStarScore.class"
printf '%s:%s\n' "$ROOT/classes" "$CP" > "$ROOT/test_classpath.txt"
printf 'base_build\t%s\npatch_build_job\t%s\nchanged_runtime_class\tKStarScore\n' "$BASE" "$SLURM_JOB_ID" > "$ROOT/build_manifest.tsv"
bash -n "$ROOT/design.slurm" "$ROOT/source/slurm/h200/node_log10_reruns_20260927.slurm"
cd "$ROOT"
find source classes -type f -print0 | sort -z | xargs -0 sha256sum > source.sha256
sha256sum cases.tsv design.slurm test_classpath.txt > launch.sha256
touch READY
echo "READY=$ROOT"
