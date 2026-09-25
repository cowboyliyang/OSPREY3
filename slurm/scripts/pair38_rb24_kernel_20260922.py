"""Slurm-only resource-aware pair ablation of the frozen RB=2/4 baseline38 runs."""
import argparse
import ast
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path('/home/users/lz280/OSPREY3-fresh-packstar')
BUILD = Path('/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588')
BASE_BUILD = Path('/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056')
BASELINE = Path('/usr/xtmp/lz280/packstar_rb38_20260918/launch_12637884/rb1')
AUDIT = Path('/usr/xtmp/lz280/packstar_pair38_20260921/audit_12677832')
JAVA_ROOT = Path('/home/users/lz280/java/jdk-17.0.2+8')
RESOURCE_FIELDS = ('cpus', 'gpus', 'mem_gib', 'heap_gib', 'host_gib')
ALLOWED_PROPERTY_CHANGES = {
    'java.io.tmpdir', 'branchdp.dp.mmap.dir', 'branchdp.dp.mmap.thresholdBytes',
    'branchdp.rootSplit.hostBudgetBytes', 'branchdp.rootSplit.gpuBudgetBytes',
    'branchdp.dp.gpu.maxGpus', 'branchdp.dp.parallel.threads',
    'packstar.dp.parallel.threads', 'packstar.pac.sampling.threads',
    'packstar.pac.frequencySeverity.tripleEtaMaxHostBytes',
    'packstar.pac.frequencySeverity.tripleEtaMaxTotalTableBytes',
    'packstar.pac.frequencySeverity.tripleEta',
    'packstar.pac.frequencySeverity.outputDir', 'osprey.bench.numCPUs',
    'osprey.bench.outputDir', 'osprey.bench.pdbPath',
}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rows(path, delimiter=','):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def keys(path):
    return {r['key']: r['value'] for r in rows(path, '\t')}


def write_table(path, data):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(data)


def checked(command, log, **kwargs):
    with log.open('w') as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True, **kwargs)


def run(root, index):
    plan = json.loads((root / 'plan.json').read_text())
    arm = plan['arm']
    assert arm in ('pair-only', 'budget-forward')
    triple_enabled = 'true' if arm == 'budget-forward' else 'false'
    case = next(c for c in plan['cases'] if c['index'] == index)
    assert int(os.environ['SLURM_CPUS_PER_TASK']) == case['cpus']
    assert os.environ.get('SLURM_JOB_ACCOUNT', 'grisman') == 'grisman'
    devices = os.environ.get('CUDA_VISIBLE_DEVICES', '').split(',')
    assert len(devices) == case['gpus'] and all(devices), devices
    gpu_text = subprocess.check_output(['nvidia-smi', '-i', ','.join(devices),
        '--query-gpu=memory.total', '--format=csv,noheader,nounits'], text=True)
    gpu_mib = [int(value.strip()) for value in gpu_text.splitlines()]
    assert len(gpu_mib) == case['gpus']
    gpu_gib = min(gpu_mib) * 80 // 100 // 1024
    env = dict(os.environ, BUILD_ROOT=plan['build'], INPUT_ROOT=str(root / 'inputs'),
               RESULT_ROOT=str(root / 'runs'), JAVA_HOME=str(JAVA_ROOT), PACKSTAR_EXPECT_GPU='')
    command = [sys.executable, str(root / 'source/portable_run.py'),
               '--cohort', 'baseline38', '--index', str(index), '--mode', 'full',
               '--arm', arm, '--seed', '42', '--rb', str(plan['rb']),
               '--gpus', str(case['gpus']), '--heap-gib', str(case['heap_gib']),
               '--host-gib', str(case['host_gib']), '--gpu-gib', str(gpu_gib)]
    task = dict(index=index, design=case['design'], job=os.environ['SLURM_JOB_ID'],
                status='RUNNING', command=command, resource_plan=case)
    task_path = root / 'tasks' / (str(index) + '.json')
    dump(task_path, task)
    proc = subprocess.run(command, env=env)
    candidates = list((root / 'runs').glob(f'{case["design"]}_full_{arm}_s42_J{os.environ["SLURM_JOB_ID"]}_T{index}'))
    assert len(candidates) == 1, candidates
    dest = candidates[0]
    manifest = json.loads((dest / 'manifest.json').read_text())
    task.update(run=str(dest), runner_exit=proc.returncode,
                java_exit=manifest.get('exit_code'), manifest_status=manifest['status'])
    if proc.returncode not in (0, 5) or manifest.get('exit_code') != 0:
        task.update(status='PROCESS_FAILED', technical_ok=False)
        dump(task_path, task)
        raise SystemExit(proc.returncode or 1)
    try:
        prior = json.loads((root / 'baseline' / case['design'] / 'manifest.json').read_text())
        assert manifest['pdb_sha256'] == prior['pdb_sha256'] == case['pdb_sha256']
        props, old_props = manifest['properties'], prior['properties']
        differences = {key: {'baseline':old_props.get(key), 'pair':props.get(key)}
            for key in set(old_props) | set(props) if old_props.get(key) != props.get(key)}
        assert not (set(differences) - ALLOWED_PROPERTY_CHANGES), differences
        assert props['branchdp.cutoff.residualBudget'] == str(plan['rb'])
        assert props['packstar.pac.frequencySeverity.proposalLearning'] == 'true'
        assert props['packstar.pac.frequencySeverity.jointMomentLearning'] == 'true'
        assert props['packstar.pac.frequencySeverity.tripleEta'] == triple_enabled
        result = rows(dest / (case['design'] + '_packstar.csv'))
        assert len(result) == case['expected_rows']
        assert {r['sequence'].strip() for r in result} == set(case['expected_sequences'])
        states = list((dest / 'adaptive_frequency_severity').glob('state-*'))
        protocols = 0
        finals = 0
        for state in states:
            protocol = keys(state / 'frequency_severity_protocol.tsv')
            if protocol:
                assert protocol['tripleEtaEnabled'] == triple_enabled
                assert protocol['proposalLearningEnabled'] == 'true'
                protocols += 1
            for path in state.glob('*selection.tsv'):
                selection = keys(path)
                if arm == 'pair-only':
                    assert int(selection.get('selectedTripleEtaPositionTriples', '0')) == 0
            final = keys(state / 'frequency_severity_final_interval.tsv')
            if final:
                if arm == 'pair-only':
                    assert int(final['selectedTripleEtaPositionTriples']) == 0
                finals += 1
        assert protocols > 0
        for row in result:
            for state in ('prot', 'lig', 'comp'):
                assert row[state + '_status'] in ('Estimated', 'Aborted', 'Unstable', 'OutOfConformations', 'OutOfLowEnergies'), row
                if row[state + '_status'] == 'Estimated':
                    assert math.isfinite(float(row[state + '_eps']))
                    assert float(row[state + '_eps']) <= 0.683001
        task.update(status='COMPLETE', technical_ok=True, rows=len(result),
                    estimated=sum(all(r[s + '_status'] == 'Estimated' for s in ('prot','lig','comp')) for r in result),
                    protocols_audited=protocols, finals_audited=finals,
                    property_differences=differences)
        dump(dest / ('pair_only_audit.json' if arm == 'pair-only' else 'optional_triple_audit.json'), task)
        dump(task_path, task)
        # Only this run's disposable DP cache; preserve all sampled/scientific evidence.
        shutil.rmtree(dest / 'dp_mmap')
        print('AUDIT_OK ' + case['design'], flush=True)
    except Exception as exc:
        task.update(status='AUDIT_FAILED', technical_ok=False, error=repr(exc))
        dump(task_path, task)
        raise


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def summarize(root):
    plan = json.loads((root / 'plan.json').read_text())
    sequences, pfuncs, systems = [], [], []
    for case in sorted(plan['cases'], key=lambda c:c['index']):
        task_file = root / 'tasks' / (str(case['index']) + '.json')
        task = json.loads(task_file.read_text()) if task_file.exists() else {}
        dest = Path(task['run']) if task.get('run') else None
        prior = root / 'baseline' / case['design']
        old = {r['sequence'].strip():r for r in rows(prior / (case['design'] + '_packstar.csv'))}
        new = {r['sequence'].strip():r for r in rows(dest / (case['design'] + '_packstar.csv'))} if dest else {}
        old_manifest = json.loads((prior / 'manifest.json').read_text())
        sysrow = dict(design=case['design'], expected=case['expected_rows'],
                      baseline_status=old_manifest['status'], pair_status=task.get('status', 'MISSING'),
                      technical_ok=task.get('technical_ok', False),
                      baseline_available=len(old), pair_available=len(new),
                      baseline_estimated=0, pair_estimated=0,
                      pair_only_success=0, triple_only_success=0,
                      node=case['node'], cpus=case['cpus'], gpus=case['gpus'],
                      memory_gib=case['mem_gib'], run=str(dest) if dest else '')
        for seq in case['expected_sequences']:
            b, p = old.get(seq), new.get(seq)
            def success(row):
                return row is not None and all(row[s+'_status']=='Estimated' for s in ('prot','lig','comp'))
            bstatus = ('Estimated' if success(b) else 'NonSuccess') if b else 'MISSING'
            pstatus = ('Estimated' if success(p) else 'NonSuccess') if p else 'MISSING'
            if not task.get('technical_ok'):
                pstatus = task.get('status', 'MISSING')
            sysrow['baseline_estimated'] += bstatus == 'Estimated'
            sysrow['pair_estimated'] += pstatus == 'Estimated'
            sysrow['pair_only_success'] += bstatus == 'NonSuccess' and pstatus == 'Estimated'
            sysrow['triple_only_success'] += bstatus == 'Estimated' and pstatus == 'NonSuccess'
            sequences.append(dict(design=case['design'], sequence=seq,
                triple_baseline=bstatus, pair_only=pstatus,
                triple_score_log10=b.get('score_log10') if b else None,
                pair_score_log10=p.get('score_log10') if p else None))
            for state in ('prot', 'lig', 'comp'):
                be = number(b.get(state + '_eps')) if b else None
                pe = number(p.get(state + '_eps')) if p else None
                bounds = [number(r.get(state + '_qstar_' + bound + '_log10')) if r else None
                          for r in (b, p) for bound in ('lb','ub')]
                both = bool(b and p and b[state+'_status']==p[state+'_status']=='Estimated' and task.get('technical_ok'))
                overlap = max(bounds[0],bounds[2]) <= min(bounds[1],bounds[3]) if both and all(v is not None for v in bounds) else None
                pfuncs.append(dict(design=case['design'], sequence=seq, state=state,
                    triple_status=b.get(state+'_status','MISSING') if b else 'MISSING',
                    pair_status=p.get(state+'_status','MISSING') if p and task.get('technical_ok') else task.get('status','MISSING'),
                    triple_epsilon=be, pair_epsilon=pe,
                    epsilon_pair_minus_triple=pe-be if both and pe is not None and be is not None else None,
                    triple_lower_log10=bounds[0], triple_upper_log10=bounds[1],
                    pair_lower_log10=bounds[2], pair_upper_log10=bounds[3], intervals_overlap=overlap))
        systems.append(sysrow)
    report = root / 'comparison'
    report.mkdir(exist_ok=True)
    write_table(report / 'systems.tsv', systems)
    write_table(report / 'sequences.tsv', sequences)
    write_table(report / 'pfuncs.tsv', pfuncs)
    transitions = {}
    for seq in sequences:
        key = seq['triple_baseline'] + ' -> ' + seq['pair_only']
        transitions[key] = transitions.get(key, 0) + 1
    summary = dict(systems=len(systems), sequences=len(sequences),
        technically_complete_systems=sum(r['technical_ok'] for r in systems),
        baseline_estimated=sum(r['baseline_estimated'] for r in systems),
        pair_estimated=sum(r['pair_estimated'] for r in systems), transitions=transitions,
        baseline_missing_sequences=sum(r['triple_baseline']=='MISSING' for r in sequences),
        no_speedup_comparison=True, rb=plan['rb'], baseline=plan['baseline'], build=plan['build'])
    dump(report / 'summary.json', summary)
    lines = [f"# RB={plan['rb']} pair-only versus prior optional-triple baseline", '',
        'Same production/benchmark source, 38 systems, actual WT + single-mutant sequence sets, seed42 and statistical settings. Pair-only disables triples. CPU/GPU counts and memory/root budgets vary by workload as requested; timing is not used to claim speedup. Fresh energy matrices in both runs.', '',
        '```json', json.dumps(summary, indent=2), '```', '',
        '| System | Baseline Estimated | Pair Estimated | Expected | Pair-only successes | Triple-only successes | Technical status |',
        '|---|---:|---:|---:|---:|---:|---|']
    for r in systems:
        lines.append('| ' + ' | '.join(str(r[k]) for k in ['design','baseline_estimated','pair_estimated','expected','pair_only_success','triple_only_success','pair_status']) + ' |')
    lines += ['', 'Only observed NonSuccess/Estimated transitions (including raw Aborted and Unstable outcomes) count as statistical recovery or loss. Missing/process-failed baseline rows are separate. Missing or interrupted baseline calculations are reported separately from statistical failures. Actual enumeration gives 21 sequences for 4wem, correcting the old runner metadata formula of 20.', '',
              'Per-sequence and per-state transitions, epsilon and interval overlap are in sequences.tsv and pfuncs.tsv. State rows can reuse the same partition function across sequences. These are paired historical outcomes with different hardware/root budgets, not repeated-seed estimates of an isolated triple effect.']
    (report / 'comparison.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary, indent=2), flush=True)
