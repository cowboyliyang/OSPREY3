"""Collect PACK* measurements without gating the independent MARK* tasks."""
import csv
import json
import os
from pathlib import Path
import re

assert os.environ.get('SLURM_JOB_ID')
root = Path(os.environ['PACK_RUN_ROOT'])
rows = list(csv.DictReader((root/'status.tsv').open(), delimiter='\t')) if (root/'status.tsv').exists() else []
results = []
for row in rows:
    output = Path(row['output']) if row['output'] else None
    wall = output/'wall.time' if output else None
    match = re.search(r'elapsed=([\d.]+)', wall.read_text()) if wall and wall.exists() else None
    results.append(dict(design=row['design'], status=row['status'],
        seconds=float(match[1]) if match else None, result_rows=row['result_rows'],
        estimated_rows=row['estimated_rows'], expected_rows=row['expected_rows'], output=row['output']))
summary = dict(measurements=results, controls_markstar_launch=False,
    note='Report speedups only from matched sequence outcomes. CPU-cohort lower tiers for 2p4a and 3bu8 use different hardware from these new GPU measurements.')
if root.exists():
    (root/'measurement_audit.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
