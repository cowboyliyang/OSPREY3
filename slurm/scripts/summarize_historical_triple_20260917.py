"""Compare complete 2xgy runs, explicitly retaining failed sequence states."""
import csv
import json
import os
import re
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')

job = sys.argv[1]
base = Path('/usr/xtmp/lz280')
root = base / f'packstar_historical_triple_20260917/A{job}'
runs = {
    'historical_four_gpu': base/'packstar_adaptive_frequency_severity/newproposal_joint_moment_dedup_38_A12509108/runs/2xgy_A12509108_T12',
    'historical_fullnode': root/f'runs/2xgy_A{job}_T0',
    'budget_forward_fullnode': base/'packstar_frontier_fullnode_20260917/A12626683/runs/2xgy_A12626683_T0',
}
results = {}
for label, path in runs.items():
    manifest = dict(line.split('\t', 1) for line in (path/'run_manifest.tsv').read_text().splitlines())
    with (path/'2xgy_packstar.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    states = ['prot_status', 'lig_status', 'comp_status']
    estimated = [row['sequence'] for row in rows if all(row[state] == 'Estimated' for state in states)]
    incomplete = [{k: row[k] for k in ['sequence'] + states}
                  for row in rows if not all(row[state] == 'Estimated' for state in states)]
    wall = float(re.search(r'elapsed=([\d.]+)', (path/'wall.time').read_text())[1])
    results[label] = dict(path=str(path), wall_s=wall, reported_s=float(rows[0]['total_time_s']),
                         rows=len(rows), all_estimated=len(estimated), estimated_sequences=sorted(estimated),
                         incomplete=incomplete, manifest=manifest)
new = results['budget_forward_fullnode']
old = results['historical_fullnode']
results['comparison'] = dict(
    current_over_historical_fullnode_wall_ratio=new['wall_s']/old['wall_s'],
    historical_resource_scale_wall_ratio=results['historical_four_gpu']['wall_s']/old['wall_s'],
    same_complete_sequence_set=set(new['estimated_sequences']) == set(old['estimated_sequences']),
    note='Whole frozen versions are compared. Different aborted-state sets make wall ratios descriptive, not matched successful-workload speedups.')
(root/'comparison.json').write_text(json.dumps(results, indent=2) + '\n')
lines = ['2xgy: complete 39-sequence P/L/PL workload', '',
         '| Version | Process seconds | Minutes | All states Estimated |',
         '|---|---:|---:|---:|']
for label in runs:
    r = results[label]
    lines.append(f"| {label} | {r['wall_s']:.2f} | {r['wall_s']/60:.2f} | {r['all_estimated']}/{r['rows']} |")
lines += ['', results['comparison']['note'], '',
          'current / old full-node wall ratio: ' + str(results['comparison']['current_over_historical_fullnode_wall_ratio'])]
(root/'comparison.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
