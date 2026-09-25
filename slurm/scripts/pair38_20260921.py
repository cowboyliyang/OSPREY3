"""Slurm-only resource-aware pair ablation of the frozen RB=1 baseline38 run."""
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


def launch():
    root = Path('/usr/xtmp/lz280/packstar_pair38_20260921') / ('launch_' + os.environ['SLURM_JOB_ID'])
    source = root / 'source'
    source.mkdir(parents=True, exist_ok=False)
    for name in ('pair38_20260921.py', 'pair38_resources_20260921.csv',
                 'run_pair38_20260921.slurm', 'summarize_pair38_20260921.slurm',
                 'launch_pair38_20260921.slurm'):
        shutil.copy2(REPO / 'slurm/scripts' / name, source / name)
    shutil.copy2(REPO / 'slurm/h200/run.py', source / 'portable_run.py')
    for path in source.glob('*.py'):
        ast.parse(path.read_text())
    for path in source.glob('*.slurm'):
        subprocess.run(['bash', '-n', str(path)], check=True)
    assert (BUILD / 'READY').is_file()
    # Current production algorithm and benchmark entry point must match baseline.
    for i, relative in enumerate(('src/main', 'src/test/java/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.java',
                                  'slurm/h200/production.properties', 'slurm/h200/baseline38.csv')):
        checked(['diff', '-qr', str(BASE_BUILD / 'source' / relative), str(BUILD / 'source' / relative)],
                root / f'baseline_source_identity_{i}.txt')
        checked(['diff', '-qr', str(REPO / relative), str(BUILD / 'source' / relative)],
                root / f'current_source_identity_{i}.txt')
    with (BUILD / 'source/slurm/h200/baseline38.csv').open() as stream:
        cohort = list(csv.reader(line for line in stream if not line.startswith('#')))
    history = {r['design']: r for r in json.loads((AUDIT / 'baseline.json').read_text())}
    allocations = rows(source / 'pair38_resources_20260921.csv')
    assert len(cohort) == len(allocations) == len(history) == 38
    assert {r[0] for r in cohort} == {r['design'] for r in allocations} == set(history)
    assert len({r['node'] for r in allocations}) == 38
    index_by_design = {r[0]: i for i, r in enumerate(cohort)}
    (root / 'inputs').mkdir()
    (root / 'preflight').mkdir()
    (root / 'runs').mkdir()
    (root / 'tasks').mkdir()
    baseline_evidence = root / 'baseline'
    baseline_evidence.mkdir()
    cases = []
    for resource in allocations:
        resource = dict(resource)
        resource.update({k: int(resource[k]) for k in RESOURCE_FIELDS})
        assert 0 < resource['host_gib'] < resource['heap_gib'] < resource['mem_gib']
        design = resource['design']
        index = index_by_design[design]
        row = cohort[index]
        prior = Path(history[design]['path'])
        prior_manifest = json.loads((prior / 'manifest.json').read_text())
        assert prior_manifest['options']['index'] == index
        assert prior_manifest['design']['mutable'] == row[5]
        assert prior_manifest['design']['flexible'] == row[6]
        pdb = Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped') / row[1] / (row[1] + '.min.reduce.renum.pdb')
        target = root / 'inputs' / row[1] / pdb.name
        target.parent.mkdir()
        shutil.copy2(pdb, target)
        digest = sha(target)
        assert digest == sha(pdb) == prior_manifest['pdb_sha256']
        evidence = baseline_evidence / design
        evidence.mkdir()
        for name in ('manifest.json', 'wall.time', 'gpus.csv', design + '_packstar.csv'):
            if (prior / name).exists():
                shutil.copy2(prior / name, evidence / name)
        cases.append(dict(resource, index=index, pdb=str(target), pdb_sha256=digest,
                          mutable=row[5], flexible=row[6], baseline=str(prior),
                          historical_elapsed_s=history[design]['elapsed_s'],
                          historical_peak_rss_gib=history[design]['peak_rss_gib'],
                          historical_status=history[design]['status']))

    def preflight(case):
        dest = root / 'preflight' / case['design']
        dest.mkdir()
        command = [str(JAVA_ROOT / 'bin/java'), '--add-modules', 'jdk.incubator.foreign',
                   '--add-opens', 'java.base/java.util=ALL-UNNAMED',
                   '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
                   '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED',
                   '-Xmx4g', '-XX:ActiveProcessorCount=2', '-Djava.io.tmpdir=' + str(dest),
                   '-Dosprey.bench.method=sequence_dump', '-Dosprey.bench.numCPUs=1',
                   '-Dosprey.bench.designId=' + case['design'],
                   '-Dosprey.bench.pdbPath=' + case['pdb'],
                   '-Dosprey.bench.outputDir=' + str(dest),
                   '-Dosprey.bench.mutable=' + case['mutable'],
                   '-Dosprey.bench.flexible=' + case['flexible'],
                   '-Dosprey.sequenceDump.maxMut=1',
                   '-Dosprey.sequenceDump.output=' + str(dest / 'sequences.tsv'),
                   '-cp', (BUILD / 'test_classpath.txt').read_text().strip(),
                   'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
        checked(command, dest / 'sequence_dump.log', cwd=dest)
        sequence_rows = rows(dest / 'sequences.tsv', '\t')
        sequences = [r['sequence'].strip() for r in sequence_rows]
        assert len(sequences) == len(set(sequences)) > 0
        prior = rows(baseline_evidence / case['design'] / (case['design'] + '_packstar.csv'))
        if prior:
            assert set(sequences) == {r['sequence'].strip() for r in prior}, case['design']
        case['expected_sequences'] = sequences
        case['expected_rows'] = len(sequences)
        print(f'PREFLIGHT {case["design"]}: {len(sequences)} sequences', flush=True)
        return case

    with ThreadPoolExecutor(max_workers=2) as pool:
        cases = list(pool.map(preflight, cases))
    baseline_systems = []
    for case in cases:
        results = rows(baseline_evidence / case['design'] / (case['design'] + '_packstar.csv'))
        estimated = sum(all(r[s + '_status'] == 'Estimated' for s in ('prot','lig','comp')) for r in results)
        baseline_systems.append(dict(design=case['design'], expected=case['expected_rows'],
            available=len(results), estimated=estimated,
            all_estimated=len(results) == case['expected_rows'] == estimated,
            pfunc_statuses={state:{status:sum(r[state + '_status']==status for r in results)
                for status in sorted({r[state + '_status'] for r in results})}
                for state in ('prot','lig','comp')}))
    baseline_summary = dict(systems=38,
        complete_systems=sum(c['all_estimated'] for c in baseline_systems),
        planned_sequences=sum(c['expected'] for c in baseline_systems),
        available_sequences=sum(c['available'] for c in baseline_systems),
        estimated_sequences=sum(c['estimated'] for c in baseline_systems),
        systems_detail=baseline_systems)
    dump(root / 'baseline_summary.json', baseline_summary)
    print('BASELINE_SUMMARY=' + json.dumps({k:v for k,v in baseline_summary.items() if k!='systems_detail'}), flush=True)
    dump(root / 'plan.json', dict(build=str(BUILD), baseline=str(BASELINE),
          rb=1, seed=42, arm='pair-only', cases=cases,
          total_sequences=sum(c['expected_rows'] for c in cases),
          scientific_comparison='completion, paired success/failure and intervals; no speedup claims',
          storage_estimate='38 private PDBs and scripts under 25 MiB; allow 100 GiB reports/EMAT, workload-dependent DP scratch up to 19 TiB and approximately 100000 files; no downloads',
          source_sha256={p.name:sha(p) for p in source.iterdir()}))
    with (root / 'scheduler_before_submission.txt').open('w') as stream:
        subprocess.run(['scontrol', 'show', 'nodes', '--oneliner'], stdout=stream, check=True)
    (root / 'READY').write_text('Current source, frozen baseline, PDB hashes and actual sequence sets validated.\n')
    jobs = []
    # Submit in descending historical duration so long jobs acquire resources first.
    with (root / 'jobs.jsonl').open('x') as registry:
        for case in sorted(cases, key=lambda c: c['historical_elapsed_s'], reverse=True):
            command = ['sbatch', '--parsable', '--account=grisman',
                       '--partition=' + case['partition'], '--nodelist=' + case['node'],
                       '--job-name=pair38_' + case['design'],
                       '--cpus-per-task=' + str(case['cpus']), '--mem=' + str(case['mem_gib']) + 'G',
                       '--gres=gpu:' + case['gpu_type'] + ':' + str(case['gpus']),
                       '--export=ALL,PAIR38_ROOT=' + str(root) + ',CASE_INDEX=' + str(case['index']),
                       str(source / 'run_pair38_20260921.slurm')]
            ident = subprocess.check_output(command, text=True).strip().split(';')[0]
            assert ident.isdigit(), ident
            record = dict(job=ident, index=case['index'], design=case['design'],
                          resources={k:case[k] for k in (*RESOURCE_FIELDS, 'node', 'partition', 'gpu_type')},
                          command=command)
            registry.write(json.dumps(record) + '\n'); registry.flush()
            jobs.append(ident)
            print('SUBMITTED ' + json.dumps(record), flush=True)
    summary = subprocess.check_output(['sbatch', '--parsable', '--account=grisman',
        '--dependency=afterany:' + ':'.join(jobs), '--export=ALL,PAIR38_ROOT=' + str(root),
        str(source / 'summarize_pair38_20260921.slurm')], text=True).strip().split(';')[0]
    dump(root / 'submission.json', dict(jobs=jobs, summary_job=summary,
         systems=len(cases), total_sequences=sum(c['expected_rows'] for c in cases),
         cpus=sum(c['cpus'] for c in cases), gpus=sum(c['gpus'] for c in cases),
         memory_gib=sum(c['mem_gib'] for c in cases)))
    (root / 'SUBMITTED').write_text('All 38 systems and afterany summary submitted.\n')
    print('ROOT=' + str(root), flush=True)
    print('SUMMARY=' + summary, flush=True)


def run(root, index):
    plan = json.loads((root / 'plan.json').read_text())
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
               '--arm', 'pair-only', '--seed', '42', '--rb', '1',
               '--gpus', str(case['gpus']), '--heap-gib', str(case['heap_gib']),
               '--host-gib', str(case['host_gib']), '--gpu-gib', str(gpu_gib)]
    task = dict(index=index, design=case['design'], job=os.environ['SLURM_JOB_ID'],
                status='RUNNING', command=command, resource_plan=case)
    task_path = root / 'tasks' / (str(index) + '.json')
    dump(task_path, task)
    proc = subprocess.run(command, env=env)
    candidates = list((root / 'runs').glob(f'{case["design"]}_full_pair-only_s42_J{os.environ["SLURM_JOB_ID"]}_T{index}'))
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
        assert props['branchdp.cutoff.residualBudget'] == '1'
        assert props['packstar.pac.frequencySeverity.proposalLearning'] == 'true'
        assert props['packstar.pac.frequencySeverity.jointMomentLearning'] == 'true'
        assert props['packstar.pac.frequencySeverity.tripleEta'] == 'false'
        result = rows(dest / (case['design'] + '_packstar.csv'))
        assert len(result) == case['expected_rows']
        assert {r['sequence'].strip() for r in result} == set(case['expected_sequences'])
        states = list((dest / 'adaptive_frequency_severity').glob('state-*'))
        protocols = 0
        finals = 0
        for state in states:
            protocol = keys(state / 'frequency_severity_protocol.tsv')
            if protocol:
                assert protocol['tripleEtaEnabled'] == 'false'
                assert protocol['proposalLearningEnabled'] == 'true'
                protocols += 1
            for path in state.glob('*selection.tsv'):
                selection = keys(path)
                assert int(selection.get('selectedTripleEtaPositionTriples', '0')) == 0
            final = keys(state / 'frequency_severity_final_interval.tsv')
            if final:
                assert int(final['selectedTripleEtaPositionTriples']) == 0
                finals += 1
        assert protocols > 0
        for row in result:
            for state in ('prot', 'lig', 'comp'):
                assert row[state + '_status'] in ('Estimated', 'Aborted'), row
                if row[state + '_status'] == 'Estimated':
                    assert math.isfinite(float(row[state + '_eps']))
                    assert float(row[state + '_eps']) <= 0.683001
        task.update(status='COMPLETE', technical_ok=True, rows=len(result),
                    estimated=sum(all(r[s + '_status'] == 'Estimated' for s in ('prot','lig','comp')) for r in result),
                    protocols_audited=protocols, finals_audited=finals,
                    property_differences=differences)
        dump(dest / 'pair_only_audit.json', task)
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
            bstatus = ('Estimated' if success(b) else 'Aborted') if b else 'MISSING'
            pstatus = ('Estimated' if success(p) else 'Aborted') if p else 'MISSING'
            if not task.get('technical_ok'):
                pstatus = task.get('status', 'MISSING')
            sysrow['baseline_estimated'] += bstatus == 'Estimated'
            sysrow['pair_estimated'] += pstatus == 'Estimated'
            sysrow['pair_only_success'] += bstatus == 'Aborted' and pstatus == 'Estimated'
            sysrow['triple_only_success'] += bstatus == 'Estimated' and pstatus == 'Aborted'
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
        no_speedup_comparison=True, baseline=str(BASELINE), build=plan['build'])
    dump(report / 'summary.json', summary)
    lines = ['# RB=1 pair-only versus prior optional-triple baseline', '',
        'Same production/benchmark source, 38 systems, actual WT + single-mutant sequence sets, seed42 and statistical settings. Pair-only disables triples. CPU/GPU counts and memory/root budgets vary by workload as requested; timing is not used to claim speedup. Fresh energy matrices in both runs.', '',
        '```json', json.dumps(summary, indent=2), '```', '',
        '| System | Baseline Estimated | Pair Estimated | Expected | Pair-only successes | Triple-only successes | Technical status |',
        '|---|---:|---:|---:|---:|---:|---|']
    for r in systems:
        lines.append('| ' + ' | '.join(str(r[k]) for k in ['design','baseline_estimated','pair_estimated','expected','pair_only_success','triple_only_success','pair_status']) + ' |')
    lines += ['', 'Only observed Aborted/Estimated transitions count as statistical recovery or loss. Missing/process-failed baseline rows are separate. In particular, baseline 3bua failed its 80-GiB host-root budget; improved availability with more memory does not establish a pair-model advantage. Actual enumeration gives 21 sequences for 4wem, correcting the old runner metadata formula of 20.', '',
              'Per-sequence and per-state transitions, epsilon and interval overlap are in sequences.tsv and pfuncs.tsv. State rows can reuse the same partition function across sequences. These are paired historical outcomes with different hardware/root budgets, not repeated-seed estimates of an isolated triple effect.']
    (report / 'comparison.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('All modes require Slurm account=grisman.')
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['launch','run','summarize'])
    parser.add_argument('--root', type=Path)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if args.root:
        assert str(args.root.absolute()).startswith('/usr/xtmp/lz280/packstar_pair38_20260921/')
    if args.mode == 'launch': launch()
    elif args.mode == 'run': run(args.root, args.index)
    else: summarize(args.root)
