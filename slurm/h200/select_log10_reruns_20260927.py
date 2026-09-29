"""Slurm-only audit of the active PRO6000 cohort; freeze observed output stalls."""
import csv
import json
import math
import os
import re
import sys
from pathlib import Path

assert os.environ['SLURM_JOB_ID']
root = Path(sys.argv[1])
old = Path('/usr/xtmp/lz280/packstar_pro6000_20260925')
status = {}
for job in (12701384, 12701385):
    for row in csv.DictReader((old / f'node_{job}/status.tsv').open(), delimiter='\t'):
        status[(str(job), int(row['index']))] = row
cases = []
fixtures = set()
pattern = re.compile(r'Z bounds: lower=(\S+), upper=(\S+), log10Lower=(\S+), log10Upper=(\S+)')
for path in sorted((old / 'results').glob('*_full_pair-only_s42_J1270138[45]_T*')):
    manifest = json.loads((path / 'manifest.json').read_text())
    index = manifest['options']['index']
    job = manifest['slurm_job_id']
    design = manifest['design']['design_id']
    csvpath = path / f'{design}_packstar.csv'
    rows = list(csv.DictReader(csvpath.open())) if csvpath.exists() and csvpath.stat().st_size else []
    audits = list((path / 'adaptive_frequency_severity').glob('*/frequency_severity_final_interval.tsv'))
    final_audit = max((p.stat().st_mtime for p in audits), default=None)
    gap = csvpath.stat().st_mtime - final_audit if rows and final_audit else None
    minimum = min((float(r[k]) for r in rows for k in r if k.endswith('_log10')
                   and r[k] and math.isfinite(float(r[k]))), default=math.inf)
    entries = []
    tail = ''
    # Stream logs only for a measured output delay or unfinished CSV.
    if not rows or (gap is not None and gap >= 60):
        with (path / 'run.log').open(errors='replace') as stream:
            for line in stream:
                if line.strip(): tail = line.strip()
                match = pattern.search(line)
                if match:
                    lo, hi, loglo, loghi = match.groups()
                    for value, log in ((lo, loglo), (hi, loghi)):
                        if math.isfinite(float(log)):
                            entries.append((value, log))
                            minimum = min(minimum, float(log))
    reason = ''
    if minimum < -324 and gap is not None and gap >= 60:
        reason = 'completed_csv_output_delay_ge_60s'
    elif minimum < -324 and csvpath.exists() and csvpath.stat().st_size == 0 and 'estimator finished' in tail:
        reason = 'all_estimators_finished_empty_csv'
    if not reason:
        continue
    assert manifest['options']['rb'] == 1 and manifest['options']['arm'] == 'pair-only'
    fixtures.update(entries)
    cases.append(dict(index=index, design_id=design, reason=reason, original_job=job,
                      original_status=status.get((job, index), {}).get('status', manifest['status']),
                      output_delay_s='' if gap is None else round(gap, 3),
                      min_log10=minimum, original_path=str(path)))
cases.sort(key=lambda r: r['index'])
assert {18, 19, 56, 57, 58, 59, 60, 61, 103, 104, 125} <= {r['index'] for r in cases}, cases
with (root / 'cases.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(cases[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(cases)
with (root / 'archived_z.tsv').open('w') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['decimal', 'log10'])
    writer.writerows(sorted(fixtures))
for row in cases:
    print(row['index'], row['design_id'], row['reason'], row['original_status'], row['output_delay_s'])
print('CASES', len(cases), 'ARCHIVED_BOUNDS', len(fixtures))
