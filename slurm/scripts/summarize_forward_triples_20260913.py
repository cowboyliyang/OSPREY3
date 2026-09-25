#!/usr/bin/env python3
"""Summarize the fixed three-case/three-arm protocol. Execute only through Slurm."""
import csv
import json
import os
import re
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run this analysis through Slurm.')
root = Path(sys.argv[1])

def table(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def keys(path):
    return {row['key']: row['value'] for row in table(path)}

rows = []
for arm in ('pair-only', 'current-selector', 'budget-forward'):
    for case in ('1gwc', '2p4a', '3k3q'):
        runs = list((root / arm / 'runs').glob(case + '_A*_T*'))
        if len(runs) != 1:
            rows.append(dict(arm=arm, case=case, status='MISSING', failure=True))
            continue
        run = runs[0]
        result = table(run / (case + '_packstar_pfunc.tsv'))
        result = result[-1] if result else {}
        log = (run / 'run.log').read_text(errors='replace') if (run / 'run.log').exists() else ''
        wall = (run / 'wall.time').read_text() if (run / 'wall.time').exists() else ''
        wall_match = re.search(r'elapsed=([\d.]+)', wall)
        totals = re.findall(r'\[PACK\*\] Total: (\d+) ms, (\d+) unique CCD calls for (\d+) sample records', log)
        states = list((run / 'adaptive_frequency_severity').glob('state-*'))
        final, winner, failure_artifact = {}, {}, {}
        for state in states:
            stage = table(state / 'frequency_severity_final.tsv')
            if stage:
                final = stage[-1]
            winner.update(keys(state / 'eta_selected.tsv'))
            for fail in state.glob('*failure*.tsv'):
                failure_artifact.update(keys(fail))
        status = result.get('status', 'MISSING')
        unique = int(totals[-1][1]) if totals else failure_artifact.get('totalCcd')
        records = int(totals[-1][2]) if totals else failure_artifact.get('totalSampleRecords')
        n = int(final['n']) if final else None
        row = dict(arm=arm, case=case, status=status, failure=status != 'Estimated',
                   wall_s=float(wall_match[1]) if wall_match else None,
                   pfunc_s=float(result['elapsed_s']) if result.get('elapsed_s') else None,
                   estimator_s=int(totals[-1][0])/1000 if totals else None,
                   final_n=n, successful_final_n=n if status == 'Estimated' else None,
                   unique_ccd=int(unique) if unique is not None else None,
                   draw_records=int(records) if records is not None else None,
                   selected_k=int(winner['tripleEtaSelectedPositionTriples']) if winner else None,
                   selected_candidate=winner.get('candidate'), epsilon=result.get('epsilon'),
                   failure_reason=failure_artifact.get('reason'), run=str(run))
        rows.append(row)

fields = ['arm', 'case', 'status', 'failure', 'wall_s', 'pfunc_s', 'estimator_s',
          'final_n', 'successful_final_n', 'unique_ccd', 'draw_records', 'selected_k',
          'selected_candidate', 'epsilon', 'failure_reason', 'run']
with (root / 'comparison.tsv').open('w') as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
summary = {}
for arm in ('pair-only', 'current-selector', 'budget-forward'):
    subset = [r for r in rows if r['arm'] == arm]
    summary[arm] = dict(cases=len(subset), failures=sum(r['failure'] for r in subset),
                        completed=sum(r['status'] != 'MISSING' for r in subset))
    for metric in ('wall_s', 'successful_final_n', 'unique_ccd', 'draw_records'):
        values = [r[metric] for r in subset if r.get(metric) is not None]
        summary[arm][metric + '_sum'] = sum(values)
        summary[arm][metric + '_count'] = len(values)
(root / 'comparison.json').write_text(json.dumps(dict(rows=rows, summary=summary), indent=2) + '\n')
lines = ['# Fixed three-arm PACK* comparison', '',
         'Single pfunc: Complex/seq0, seed 42. Missing/aborted final stages are not counted as sample savings.', '',
         '| Case | Arm | Status | Wall s | Final N | Unique CCD | Draws | K |',
         '|---|---|---|---:|---:|---:|---:|---:|']
for r in rows:
    lines.append('| ' + ' | '.join(str(r.get(k)) if r.get(k) is not None else 'NA'
                                  for k in ('case', 'arm', 'status', 'wall_s', 'final_n', 'unique_ccd', 'draw_records', 'selected_k')) + ' |')
lines += ['', 'All nine runs share one tested build. Each case runs all arms on one allocated node.',
          'This small comparison does not estimate performance or failure rates over the 38-design population.']
(root / 'comparison.md').write_text('\n'.join(lines) + '\n')
print(json.dumps(summary, indent=2))
