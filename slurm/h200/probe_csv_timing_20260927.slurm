#!/bin/bash
#SBATCH --job-name=ps-csv-timing
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_csv_timing_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_csv_timing_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
export ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/csv_timing_${SLURM_JOB_ID}
mkdir "$ROOT"
mkdir "$ROOT/tmp" "$ROOT/classes" "$ROOT/exports"
export TMPDIR="$ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
BASE=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/complete_12706035
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/CsvOutputTimingProbe.java "$ROOT/"
cp "$BASE/all_162.tsv" "$ROOT/all_162.tsv"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
CP=$(cat "$BASE/test_classpath.txt")
JDK=/home/users/lz280/java/jdk-17.0.2+8/bin
"$JDK/javac" -implicit:none -cp "$CP" -d "$ROOT/classes" "$ROOT/CsvOutputTimingProbe.java"
timeout --kill-after=5s 180s "$JDK/java" --add-modules jdk.incubator.foreign -Xmx2g \
    -Djava.io.tmpdir="$TMPDIR" -cp "$ROOT/classes:$CP" CsvOutputTimingProbe "$ROOT/all_162.tsv" "$ROOT/exports" \
    > "$ROOT/formatter.log" 2>&1
python3 - <<'PY'
import csv, json, os, statistics
from pathlib import Path
root = Path(os.environ['ROOT'])
audit = list(csv.DictReader((root / 'all_162.tsv').open(), delimiter='\t'))
timings = list(csv.DictReader((root / 'exports/csv_timings.tsv').open(), delimiter='\t'))
completed = [r for r in audit if r['stage'] == 'CSV_WRITTEN']
assert len(timings) == len(completed)
not_selected = [r for r in completed if r['any_extreme'] == 'False']
print('TIMING_FIXTURES: rounded exported bounds reconstructed as BigDecimals; real sorter, formatter, filesystem write and close; not a scientific rerun')
print('COMPLETED_REPLAYED', len(timings), 'ROWS', sum(int(r['rows']) for r in timings))
for name in ('first_export_ms', 'repeat_export_ms'):
    vals = [float(r[name]) for r in timings]
    print(name, 'min', min(vals), 'median', statistics.median(vals), 'max', max(vals), 'over500ms', sum(v > 500 for v in vals))
print('UNSELECTED_COMPLETED_CASES', len(not_selected), 'MAX_ORIGINAL_GAP_S', max(float(r['output_gap_s']) for r in not_selected))
for row in sorted(not_selected, key=lambda r: float(r['output_gap_s']), reverse=True)[:12]:
    print('UNSELECTED', row['index'], row['design_id'], 'original_gap_s', row['output_gap_s'], 'min_log10', row['min_log10'], 'max_log10', row['max_log10'])
print('SLOWEST_REPLAYS')
for row in sorted(timings, key=lambda r: float(r['first_export_ms']), reverse=True)[:10]: print(json.dumps(row))
PY
