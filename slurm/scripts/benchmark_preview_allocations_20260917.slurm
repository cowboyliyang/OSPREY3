#!/bin/bash
#SBATCH --job-name=ps_preview_allocations
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${NEW_BUILD:?}"
BASE_BUILD=/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626676
ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/preview_${SLURM_JOB_ID}
mkdir -p "$ROOT/classes" "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
JAVA_ROOT=/home/users/lz280/java/jdk-17.0.2+8
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/DecompositionPreviewSpeedBench.java "$ROOT/"
"$JAVA_ROOT/bin/javac" -cp "$(cat "$NEW_BUILD/test_classpath.txt")" -d "$ROOT/classes" "$ROOT/DecompositionPreviewSpeedBench.java"
for mode in base optimized; do
    build="$BASE_BUILD"
    if [[ "$mode" == optimized ]]; then build="$NEW_BUILD"; fi
    "$JAVA_ROOT/bin/java" -Xmx6g -XX:ActiveProcessorCount=1 -Djava.io.tmpdir="$TMPDIR" \
        -cp "$ROOT/classes:$(cat "$build/test_classpath.txt")" \
        edu.duke.cs.osprey.branchdp.DecompositionPreviewSpeedBench > "$ROOT/$mode.tsv"
done
python3 - "$ROOT" <<'PY'
from pathlib import Path
import json, sys
root = Path(sys.argv[1])
tables = [[line.split('\t') for line in (root/f'{mode}.tsv').read_text().splitlines()
           if line.startswith('RESULT\t')] for mode in ('base', 'optimized')]
assert len(tables[0]) == len(tables[1]) == 9
report = []
for old, new in zip(*tables):
    assert old[1] == new[1] and old[3] == new[3], (old, new)
    a, b = float(old[2]), float(new[2])
    report.append(dict(trial=int(old[1]), exact_ordered_tree=True,
                       old_s=a, new_s=b, speedup=a/b))
result = dict(cases=report, summed_speedup=sum(r['old_s'] for r in report)/sum(r['new_s'] for r in report))
(root/'comparison.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
PY
