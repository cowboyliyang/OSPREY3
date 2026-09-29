"""Include every observed value using the patched extreme-exponent output path."""
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

assert os.environ['SLURM_JOB_ID']
root = Path(sys.argv[1])
build = Path(sys.argv[2])
rows = list(csv.DictReader((root / 'all_162.tsv').open(), delimiter='\t'))
selected = [r for r in rows if r['any_extreme'] == 'True']
assert selected and {18, 19, 56, 57, 58, 59, 60, 61, 103, 104, 125} <= {int(r['index']) for r in selected}
config = build / 'source/slurm/h200'
runner_sha = hashlib.sha256((config / 'run.py').read_bytes()).hexdigest()
manifest_sha = hashlib.sha256((config / 'frontier_active.tsv').read_bytes()).hexdigest()
expected = dict(cohort='frontier', mode='full', arm='pair-only', seed=42, rb=1,
                gpus=4, heap_gib=850, host_gib=800, gpu_gib=85)
production = dict(line.split('=', 1) for line in (config / 'production.properties').read_text().splitlines() if line and not line.startswith('#'))
production.update({'branchdp.cutoff.residualBudget': '1', 'packstar.pac.frequencySeverity.tripleEta': 'false'})
for row in selected:
    m = json.loads((Path(row['path']) / 'manifest.json').read_text())
    assert all(m['options'][k] == v for k, v in expected.items()), row
    assert m['cpus'] == 128 and m['precision'] == 'FP64' and m['ccd'] == 'CPU'
    assert m['runner_sha256'] == runner_sha and m['frontier_manifest_sha256'] == manifest_sha
    assert all(m['properties'][k] == v for k, v in production.items()), row
with (root / 'cases.tsv').open('w') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(selected)
print('SELECTED_CASES', len(selected))
print('RULE: finite log10 below -308 or at least 309 in Z/K* output, regardless of observed delay')
print('PASS: runner, positions/sequences, all production properties and resource options match')
print('UNRESOLVED_WITHOUT_OBSERVATIONS', sum(r['stage'] in ('NOT_STARTED', 'COMPUTING', 'ENDED_WITHOUT_CSV') for r in rows))
for r in selected:
    print(r['index'], r['design_id'], r['stage'], 'gap_s=', r['output_gap_s'], 'min_log10=', r['min_log10'])
