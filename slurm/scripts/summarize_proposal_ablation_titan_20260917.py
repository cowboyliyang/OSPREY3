"""Summarize all predeclared runs, including missing/aborted/time-limited arms."""
import csv
import json
import os
import re
import statistics
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')
root = Path(sys.argv[1])
systems = ('1gwc', '2p4a', '2xxm', '3k3q', '2xgy', '3ma2', '4wyu', '5a6y')
arms = ('no-learning', 'pair-only', 'decomposition-cost', 'budget-forward')

def table(path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def keys(path):
    return {row['key']: row['value'] for row in table(path)}

def manifest(path):
    if not path.exists():
        return {}
    return dict(line.split('\t', 1) for line in path.read_text().splitlines() if '\t' in line)

def number(value, kind=float):
    return kind(value) if value not in (None, '', 'NA') else None

rows = []
for case in systems:
    for seed in (42, 43, 44):
        for arm in arms:
            run = root / case / f'seed_{seed}' / arm
            results = table(run / f'{case}_packstar_pfunc.tsv')
            result = results[-1] if results else {}
            meta = manifest(run / 'run_manifest.tsv')
            log = (run / 'run.log').read_text(errors='replace') if (run / 'run.log').exists() else ''
            wall = (run / 'wall.time').read_text() if (run / 'wall.time').exists() else ''
            time_match = re.search(r'elapsed=([\d.]+) maxRssKiB=(\d+)', wall)
            totals = re.findall(r'\[PACK\*\] Total: (\d+) ms, (\d+) unique CCD calls for (\d+) sample records', log)
            final, winner, failure = {}, {}, {}
            for state in (run / 'adaptive_frequency_severity').glob('state-*'):
                stage = table(state / 'frequency_severity_final.tsv')
                if stage:
                    final = stage[-1]
                winner.update(keys(state / 'eta_selected.tsv'))
                for path in state.glob('*failure*.tsv'):
                    failure.update(keys(path))
            exit_code = meta.get('javaStatus')
            status = result.get('status', 'MISSING')
            if exit_code == '124':
                status = 'TIMEOUT'
            elif exit_code not in (None, '0'):
                status = f'PROCESS_EXIT_{exit_code}'
            elif not results and run.exists():
                status = 'INCOMPLETE'
            rows.append(dict(case=case, seed=seed, arm=arm, status=status,
                node=meta.get('node'), arm_order=number(meta.get('armOrder'), int),
                java_exit=exit_code, wall_s=float(time_match[1]) if time_match else None,
                max_rss_kib=int(time_match[2]) if time_match else None,
                pfunc_s=number(result.get('elapsed_s')),
                estimator_s=int(totals[-1][0])/1000 if totals else None,
                final_n=number(final.get('n'), int),
                unique_ccd=int(totals[-1][1]) if totals else number(failure.get('totalCcd'), int),
                draws=int(totals[-1][2]) if totals else number(failure.get('totalSampleRecords'), int),
                selected_k=number(winner.get('tripleEtaSelectedPositionTriples'), int),
                corrected_dp_sweeps=number(winner.get('proposalDpSweeps'), int),
                selected_candidate=winner.get('candidate'),
                epsilon=number(result.get('epsilon')), lower_log10=number(result.get('lower_log10')),
                upper_log10=number(result.get('upper_log10')),
                failure_reason=failure.get('reason'), run=str(run)))

with (root / 'comparison.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
summary = {}
for arm in arms:
    selected = [r for r in rows if r['arm'] == arm]
    statuses = sorted({r['status'] for r in selected})
    summary[arm] = {'attempts_planned': len(selected),
                    'statuses': {s: sum(r['status'] == s for r in selected) for s in statuses}}
paired = []
for case in systems:
    for seed in (42, 43, 44):
        group = {r['arm']: r for r in rows if r['case'] == case and r['seed'] == seed}
        base = group['pair-only']
        for arm in arms:
            other = group[arm]
            if base['status'] == other['status'] == 'Estimated' and base['wall_s'] and other['wall_s']:
                paired.append(dict(case=case, seed=seed, arm=arm,
                                   wall_relative_to_pair=other['wall_s']/base['wall_s']))
(root / 'comparison.json').write_text(json.dumps(dict(rows=rows, summary=summary, matched_successful_pairs=paired), indent=2) + '\n')
lines = ['# Proposal ablation on Titan V', '',
         '8 systems x 3 seeds x 4 arms; Complex/seq0 only. Every system/seed runs all arms on one exclusive 48-CPU/four-Titan-V node. FP64; epsilon=0.683. EMAT preparation excluded.', '',
         'No-learning retains calibration/sizing/validation/monitor/final budgets and the conditional S0=20 premise, but fixes the proposal to the original DP model. Missing/failed/aborted arms remain in the denominator.', '',
         '| Arm | Estimated / 24 | Status counts | Median wall / pair-only on matched successes |',
         '|---|---:|---|---:|']
for arm in arms:
    ratios = [r['wall_relative_to_pair'] for r in paired if r['arm'] == arm]
    ratio = f'{statistics.median(ratios):.3f} (n={len(ratios)})' if ratios else 'NA'
    lines.append(f"| {arm} | {summary[arm]['statuses'].get('Estimated', 0)} | {summary[arm]['statuses']} | {ratio} |")
lines += ['', '| System | Seed | Arm | Status | Wall s | Final N | Unique CCD | K | Peak RSS KiB |',
          '|---|---:|---|---|---:|---:|---:|---:|---:|']
for row in rows:
    lines.append('| ' + ' | '.join(str(row[k]) if row[k] is not None else 'NA' for k in
        ('case', 'seed', 'arm', 'status', 'wall_s', 'final_n', 'unique_ccd', 'selected_k', 'max_rss_kib')) + ' |')
lines += ['', 'Peak RSS measures the whole process, not GPU memory. Structural DP estimates, fit diagnostics and phase logs remain in each run directory. Successful-pair time ratios must be read with completion rates; no main configuration is selected automatically.']
(root / 'comparison.md').write_text('\n'.join(lines) + '\n')
print(json.dumps(summary, indent=2))
