#!/usr/bin/env python3
"""Slurm-only summary. A row succeeds only when all three states Estimated."""
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
systems = ('2xgy', '4u3s', '1a0r', '3ma2', '4wyu', '5a6y')
summary = []
for task, system in enumerate(systems):
    run = root / 'runs' / f'{system}_{root.name}_T{task}'
    manifest_file = run / 'run_manifest.tsv'
    manifest = dict(line.split('\t', 1) for line in manifest_file.read_text().splitlines()) if manifest_file.exists() else {}
    result = run / f'{system}_packstar.csv'
    rows = list(csv.DictReader(result.open())) if result.exists() else []
    statuses = Counter(tuple(row.get(f'{state}_status', 'missing') for state in ('prot', 'lig', 'comp')) for row in rows)
    elapsed = sorted({float(row['total_time_s']) for row in rows})
    wall_file = run / 'wall.time'
    wall_match = re.search(r'elapsed=([0-9.]+)', wall_file.read_text()) if wall_file.exists() else None
    all_estimated = statuses.get(('Estimated', 'Estimated', 'Estimated'), 0)
    summary.append(dict(
        system=system, task=task, java_status=manifest.get('javaStatus', 'missing'),
        node=manifest.get('node', 'missing'), rows=len(rows), all_estimated=all_estimated,
        other_rows=len(rows)-all_estimated,
        total_time_s=elapsed[0] if len(elapsed) == 1 else None,
        process_wall_s=float(wall_match[1]) if wall_match else None,
        state_status_counts={'/'.join(key): count for key, count in statuses.items()},
        result=str(result), run_dir=str(run)))
root.mkdir(parents=True, exist_ok=True)
(root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
fields = ('system', 'task', 'java_status', 'node', 'rows', 'all_estimated', 'other_rows', 'total_time_s', 'process_wall_s')
with (root / 'summary.tsv').open('w') as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(summary)
lines = ['# PACK* FP64 full-node baseline', '',
         '104 CPU threads, 8 A5000 GPUs, all node memory; original flexibility, budget-forward.',
         'EMAT inputs reuse historical prepared caches. Times are not fresh-EMAT end-to-end timings.',
         'A row is successful only when P, L and PL are all Estimated. Missing results and aborted states remain visible.', '',
         '| System | Node | Java exit | All Estimated / rows | PACK seconds | Process seconds |',
         '|---|---|---:|---:|---:|---:|']
for row in summary:
    lines.append(f"| {row['system']} | {row['node']} | {row['java_status']} | {row['all_estimated']}/{row['rows']} | {row['total_time_s']} | {row['process_wall_s']} |")
(root / 'summary.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
