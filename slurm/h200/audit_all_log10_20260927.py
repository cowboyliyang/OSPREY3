"""Audit every design in the active PRO6000 manifest, including unstarted designs."""
import collections
import csv
import json
import math
import os
import re
import sys
import time
from pathlib import Path

assert os.environ['SLURM_JOB_ID']
root = Path(sys.argv[1])
old = Path('/usr/xtmp/lz280/packstar_pro6000_20260925')
designs = list(csv.DictReader((old / 'duke_12701382/source/slurm/h200/frontier_active.tsv').open(), delimiter='\t'))
selected = {int(r['index']) for r in csv.DictReader(Path('/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/build_12706032/cases.tsv').open(), delimiter='\t')}
history = collections.defaultdict(list)
for r in csv.DictReader(Path('/usr/xtmp/lz280/pair_unstable_audit_20260927_12706008/runs.tsv').open(), delimiter='\t'):
    if r['rb'] == '1': history[r['design']].append(r)
statuses = {}
for job in (12701384, 12701385):
    for row in csv.DictReader((old / f'node_{job}/status.tsv').open(), delimiter='\t'):
        statuses[int(row['index'])] = row
manifests = {}
for path in (old / 'results').glob('*_full_pair-only_s42_J1270138[45]_T*/manifest.json'):
    m = json.loads(path.read_text())
    manifests[m['options']['index']] = (path.parent, m)
pattern = re.compile(r'Z bounds: lower=(\S+), upper=(\S+), log10Lower=(\S+), log10Upper=(\S+)')
records = []
bound_rows = []
for i, design in enumerate(designs):
    assert i == int(design['task_id'])
    historical = history[design['system']]
    baseline_min = min((float(r['minimum_log10']) for r in historical if r['minimum_log10']), default=None)
    baseline_gap = max((float(r['output_gap_s']) for r in historical if r['output_gap_s']), default=None)
    record = dict(index=i, design_id=design['design_id'], system=design['system'],
                  stage='NOT_STARTED', original_status='', rows=0, expected=int(design['expected_sequences']),
                  output_gap_s=None, csv_present=False, min_log10=None, max_log10=None,
                  negative_extreme_sequences=0, any_extreme_sequences=0,
                  negative_extreme=False, any_extreme=False, empty_csv_after_estimation=False,
                  already_queued=i in selected, baseline_min_log10=baseline_min, baseline_gap_s=baseline_gap,
                  path='')
    if i in manifests:
        path, m = manifests[i]
        record['path'] = str(path)
        record['original_status'] = statuses.get(i, {}).get('status', m['status'])
        csvpath = path / f"{design['design_id']}_packstar.csv"
        record['csv_present'] = csvpath.exists()
        rows = list(csv.DictReader(csvpath.open())) if csvpath.exists() and csvpath.stat().st_size else []
        record['rows'] = len(rows)
        values = []
        for row in rows:
            vals = [float(v) for k, v in row.items() if k.endswith('_log10') and v and math.isfinite(float(v))]
            values.extend(vals)
            record['negative_extreme_sequences'] += any(v < -324 for v in vals)
            record['any_extreme_sequences'] += any(v < -308 or v >= 309 for v in vals)
        tail = ''
        if (path / 'run.log').exists():
            for line in (path / 'run.log').open(errors='replace'):
                if line.strip(): tail = line.strip()
                match = pattern.search(line)
                if match:
                    lo, hi, loglo, loghi = match.groups()
                    for z, log in ((lo, loglo), (hi, loghi)):
                        if math.isfinite(float(log)):
                            values.append(float(log))
                            bound_rows.append(dict(index=i, design_id=design['design_id'], decimal=z, log10=log))
        if values:
            record['min_log10'] = min(values)
            record['max_log10'] = max(values)
            record['negative_extreme'] = min(values) < -324
            record['any_extreme'] = min(values) < -308 or max(values) >= 309
        audits = list((path / 'adaptive_frequency_severity').glob('*/frequency_severity_final_interval.tsv'))
        last_audit = max((p.stat().st_mtime for p in audits), default=None)
        if rows:
            record['stage'] = 'CSV_WRITTEN'
            if last_audit:
                record['output_gap_s'] = round(csvpath.stat().st_mtime - last_audit, 3)
        elif csvpath.exists() and csvpath.stat().st_size == 0 and 'estimator finished' in tail:
            record['stage'] = 'OUTPUT_STALLED'
            record['empty_csv_after_estimation'] = True
        else:
            record['stage'] = 'COMPUTING' if i not in statuses else 'ENDED_WITHOUT_CSV'
    records.append(record)

def write(name, data):
    with (root / name).open('w') as stream:
        w = csv.DictWriter(stream, fieldnames=list(data[0]), delimiter='\t')
        w.writeheader()
        w.writerows(data)

write('all_162.tsv', records)
write('all_logged_bounds.tsv', bound_rows)
missing_slow = [r for r in records if not r['already_queued'] and r['negative_extreme'] and r['output_gap_s'] is not None and r['output_gap_s'] >= 10]
extra_extreme = [r for r in records if not r['already_queued'] and r['negative_extreme']]
pending_risk = [r for r in records if r['stage'] in ('NOT_STARTED', 'COMPUTING') and r['baseline_min_log10'] is not None and r['baseline_min_log10'] < -10000]
write('additional_extreme.tsv', extra_extreme)
if missing_slow: write('additional_delay_ge10s.tsv', missing_slow)
if pending_risk: write('unfinished_with_extreme_baseline.tsv', pending_risk)
summary = dict(audit_epoch=time.time(), total=len(records), stages=dict(collections.Counter(r['stage'] for r in records)),
               negative_extreme_cases=sum(r['negative_extreme'] for r in records),
               any_extreme_cases=sum(r['any_extreme'] for r in records),
               output_gaps={str(t): sum(r['output_gap_s'] is not None and r['output_gap_s'] >= t for r in records) for t in (1, 5, 10, 30, 60)},
               previously_queued=len(selected), additional_negative_extreme=len(extra_extreme),
               additional_delay_ge10s=len(missing_slow), unfinished_extreme_baseline=len(pending_risk))
(root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
for label, data in [('ADDITIONAL_DELAY_GE10S', missing_slow), ('UNFINISHED_WITH_EXTREME_BASELINE', pending_risk)]:
    print(label)
    for r in data:
        print(r['index'], r['design_id'], r['stage'], 'gap_s=', r['output_gap_s'], 'min_log10=', r['min_log10'], 'baseline_gap_s=', r['baseline_gap_s'])

