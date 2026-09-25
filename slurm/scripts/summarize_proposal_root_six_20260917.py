#!/usr/bin/env python3
"""Slurm-only audit of the six proposal-root reruns against the latest baseline."""
import csv
import json
import re
import sys
from pathlib import Path


def read_run(run, system):
    csv_path = run / f'{system}_packstar.csv'
    rows = list(csv.DictReader(csv_path.open())) if csv_path.exists() else []
    metadata_path = run / 'run_manifest.tsv'
    metadata = dict(line.split('\t', 1) for line in metadata_path.read_text().splitlines()) if metadata_path.exists() else {}
    wall_path = run / 'wall.time'
    elapsed = re.search(r'elapsed=([\d.]+)', wall_path.read_text()) if wall_path.exists() else None
    hashes = {}
    if (run / 'input.sha256').exists():
        for line in (run / 'input.sha256').read_text().splitlines():
            digest, path = line.split(maxsplit=1)
            hashes[Path(path).name] = digest
    roots = []
    log_path = run / 'run.log'
    if log_path.exists():
        with log_path.open(errors='replace') as log:
            for line in log:
                if '[PACK*-triple-eta-dp] rebuilt exact proposal' in line:
                    match = re.search(r'rootSplit=(\d+), rootSelectionMs=([\d.Ee+-]+), gpuWork=(\d+)', line)
                    if match:
                        roots.append(dict(split=int(match[1]), milliseconds=float(match[2]), work=int(match[3])))
    success = sorted(r['sequence'] for r in rows if all(r[f'{state}_status'] == 'Estimated' for state in ('prot', 'lig', 'comp')))
    return dict(run=str(run), wall_s=float(elapsed[1]) if elapsed else None,
                rows=len(rows), success=success, sequences=sorted(r['sequence'] for r in rows),
                metadata=metadata, input_sha256=hashes,
                proposal_root_selections=len(roots), nonzero_proposal_roots=sum(r['split'] != 0 for r in roots),
                proposal_root_seconds=sum(r['milliseconds'] for r in roots) / 1000,
                proposal_roots=roots)


def main():
    job = sys.argv[1]
    root = Path('/usr/xtmp/lz280/packstar_proposal_root_20260917') / f'A{job}'
    baseline_root = Path('/usr/xtmp/lz280/packstar_frontier_fullnode_20260917')
    systems = ['2xgy', '4u3s', '1a0r', '3ma2', '4wyu', '5a6y']
    comparisons = []
    resource_keys = ['cpuThreads', 'heapGiB', 'hostBudgetGiB', 'gpuName', 'branchDpGpuCount',
                     'gpuBudgetBytes', 'precision', 'ccdMode', 'targetEpsilon', 'flexAdded']
    for task, system in enumerate(systems):
        baseline_job = '12626765' if task == 0 else '12627258'
        old = read_run(baseline_root / f'A{baseline_job}/runs/{system}_A{baseline_job}_T{task}', system)
        new = read_run(root / f'runs/{system}_A{job}_T{task}', system)
        comparisons.append(dict(system=system, task=task, baseline_job=baseline_job,
            baseline=old, rerun=new,
            same_sequence_workload=bool(new['sequences']) and new['sequences'] == old['sequences'],
            same_inputs=bool(new['input_sha256']) and new['input_sha256'] == old['input_sha256'],
            resource_differences={k: [old['metadata'].get(k), new['metadata'].get(k)] for k in resource_keys
                                  if old['metadata'].get(k) != new['metadata'].get(k)},
            lost_successes=sorted(set(old['success']) - set(new['success'])),
            gained_successes=sorted(set(new['success']) - set(old['success'])),
            wall_ratio=old['wall_s'] / new['wall_s'] if old['wall_s'] and new['wall_s'] else None))
    root.mkdir(parents=True, exist_ok=True)
    (root / 'six_system_comparison.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    lines = ['Proposal-root optimization: six systems, 8 A5000 / 104 CPU per task.', '',
             '| System | Before (s) | After (s) | Before/after | Estimated before | Estimated after | Nonzero proposal roots |',
             '|---|---:|---:|---:|---:|---:|---:|']
    def fmt(x):
        return f'{x:.2f}' if x is not None else 'missing'
    for c in comparisons:
        a, b = c['baseline'], c['rerun']
        lines.append(f"| {c['system']} | {fmt(a['wall_s'])} | {fmt(b['wall_s'])} | {fmt(c['wall_ratio'])} | "
                     f"{len(a['success'])}/{a['rows']} | {len(b['success'])}/{b['rows']} | "
                     f"{b['nonzero_proposal_roots']}/{b['proposal_root_selections']} |")
    lines += ['', 'Whole-run timings include learning, root selection and DP; EMAT is copied from the same historical cache.',
              'Root changes can alter selected proposals, samples and success sets. Inspect the JSON before interpreting timing ratios.']
    for c in comparisons:
        complete = c['rerun']['metadata'].get('javaStatus') == '0' and c['same_sequence_workload']
        lines.append(f"{c['system']}: complete={complete}, same_inputs={c['same_inputs']}, "
                     f"resource_differences={c['resource_differences']}, "
                     f"lost_successes={c['lost_successes']}, gained_successes={c['gained_successes']}")
    (root / 'six_system_comparison.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    if any(c['rerun']['metadata'].get('javaStatus') != '0' or not c['same_sequence_workload']
           or not c['same_inputs'] or c['resource_differences'] for c in comparisons):
        raise SystemExit('Incomplete or mismatched runs: inspect six_system_comparison.json')


if __name__ == '__main__':
    main()
