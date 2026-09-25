"""Prepare, execute, audit and summarize the predeclared matched experiment.

All entry points require Slurm. No algorithm implementation is changed here.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys

JAVA = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'
BENCH = 'edu.duke.cs.osprey.markstar.bench.GenericPDBBench'
SPEC = Path(__file__).with_name('hardcase_cases_20260921.json')
OPEN = ['--add-opens', 'java.base/java.util=ALL-UNNAMED', '--add-opens',
        'java.base/java.lang=ALL-UNNAMED', '--add-opens',
        'java.base/java.lang.invoke=ALL-UNNAMED']


def table(path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def keys(path):
    return {r['key']: r['value'] for r in table(path)}


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2) + '\n')


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def checked_run(command, log, **kwargs):
    with log.open('w') as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, **kwargs)


def java_command(build, props, heap='280g'):
    return [JAVA, *OPEN, '-Xmx' + heap, '-Xms1g', '-XX:-UseSuperWord',
            '-XX:ActiveProcessorCount=48',
            *['-D' + k + '=' + str(v) for k, v in props.items()],
            '-cp', (build / 'test_classpath.txt').read_text().strip(), BENCH]


def prepare(build):
    spec = json.loads(SPEC.read_text())
    inputs = build / 'inputs'
    inputs.mkdir()
    csv_path = Path('/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv')
    rows = {}
    for row in csv.reader(csv_path.read_text().splitlines()):
        if row and row[0] in {c['design'] for c in spec['cases']}:
            rows[row[0]] = row
    input_hashes = {}
    sequence_maps = {}
    for design, row in rows.items():
        dest = inputs / design
        dest.mkdir()
        pdb = Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped') / row[1] / (row[1] + '.min.reduce.renum.pdb')
        sources = {dest / 'input.pdb': pdb}
        for state in ('complex', 'protein', 'ligand'):
            for old, new in (('rigid', 'rigid'), ('min', 'minimizing')):
                source = Path('/usr/xtmp/lz280/bench_comparison/results/emat_cache') / design / f'export.{state}.{old}.dat'
                sources[dest / f'packstar.{state}.{new}.dat'] = source
        for target, source in sources.items():
            shutil.copy2(source, target)
            digest = sha(source)
            assert sha(target) == digest, source
            input_hashes[str(target.relative_to(inputs))] = {'source': str(source), 'sha256': digest}
        (dest / 'design_row.csv').write_text(','.join(row) + '\n')
        props = {'osprey.bench.method': 'sequence_dump', 'osprey.bench.outputDir': dest,
                 'osprey.bench.designId': design, 'osprey.bench.pdbPath': dest / 'input.pdb',
                 'osprey.bench.mutable': row[5], 'osprey.bench.flexible': row[6],
                 'osprey.bench.numCPUs': 1, 'osprey.sequenceDump.maxMut': 1,
                 'osprey.sequenceDump.output': dest / 'sequences.tsv',
                 'java.io.tmpdir': build / 'tmp'}
        checked_run(java_command(build, props, '8g'), dest / 'sequence_dump.log', cwd=build / 'source')
        sequence_maps[design] = {int(r['seq_index']): r['sequence'].strip() for r in table(dest / 'sequences.tsv')}
    assert len(rows) == len({c['design'] for c in spec['cases']})
    old_cases_path = Path('/usr/xtmp/lz280/packstar_rescue_20260905/selection_12505406/candidates.tsv')
    rescued_path = Path('/usr/xtmp/lz280/packstar_hardcase_writing_audit_20260920/12659842/rescue_metadata.tsv')
    old_cases = {r['case_id']: r for r in table(old_cases_path)}
    rescued = {r['case_id']: r for r in table(rescued_path)}
    historical = []
    for case in spec['cases']:
        assert sequence_maps[case['design']][case['seq_index']] == case['sequence'], case
        case['mutable'], case['flexible'] = rows[case['design']][5:7]
        old, new = old_cases[case['rescue_case_id']], rescued[case['rescue_case_id']]
        assert old['packstar_status'] == 'Aborted', case
        for key in ('design', 'state', 'sequence'):
            assert old[key] == new[key] == case[key], (case, key)
        assert int(new['sequence_index']) == case['seq_index']
        assert old['mutable'] == case['mutable'] and old['flexible'] == case['flexible']
        final_path = Path(case['historical_final_path'])
        assert final_path.parent.name == Path(new['state_path']).name, case
        final = keys(final_path)
        assert final['selectedTripleEtaActive'] == 'true' and int(final['selectedTripleEtaPositionTriples']) > 0
        assert final['targetReached'] == 'true' and final['certificateValid'] == 'true'
        prior_run = final_path.parents[2]
        prior_single = table(prior_run / 'result.tsv')
        if prior_single:
            matches = [r for r in prior_single if r['sequence'] == case['sequence'] and r['state'].lower() == case['state']]
            assert len(matches) == 1 and matches[0]['status'] == 'Estimated'
        else:
            with (prior_run / (case['design'] + '_packstar.csv')).open() as stream:
                matches = [r for r in csv.DictReader(stream) if r['sequence'] == case['sequence']]
            assert len(matches) == 1 and matches[0][{'complex': 'comp', 'protein': 'prot', 'ligand': 'lig'}[case['state']] + '_status'] == 'Estimated'
        historical.append({'case': case['id'], 'old': old, 'historical_final_path': str(final_path),
                           'historical_final_sha256': sha(final_path), 'final': final, 'result': matches[0]})
    dump(build / 'cases.json', spec)
    dump(build / 'input_manifest.json', input_hashes)
    dump(build / 'historical_evidence.json', {'old_manifest': str(old_cases_path),
         'old_manifest_sha256': sha(old_cases_path), 'rescue_metadata': str(rescued_path),
         'rescue_metadata_sha256': sha(rescued_path), 'cases': historical})
    print(f'PREPARED: current enumeration and historical final-triple success verified for {len(spec["cases"])} cases', flush=True)


def common_properties(case, seed, run):
    fixed = '''
branchdp.dp.cache=false
packstar.dp.cache=false
branchdp.dp.tableMode=auto_mmap
branchdp.dp.mmap.thresholdBytes=137438953472
branchdp.dp.mmap.skipInitialFill=true
branchdp.decomp.strategy=weighted_hicks
branchdp.decomp.weightedHicks.restarts=24
branchdp.decomp.weightedHicks.randomMoves=32
branchdp.dp.parallel.threads=48
packstar.dp.parallel.threads=48
branchdp.dp.gpu=true
branchdp.dp.gpu.failIfNoGpuPath=true
branchdp.dp.gpu.failIfExceedsVram=true
branchdp.dp.gpu.multiGpu=true
branchdp.dp.gpu.maxGpus=4
branchdp.dp.gpu.minMStatesPerGpu=4096
branchdp.dp.gpu.childSliceMaxBytes=2147483648
branchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296
branchdp.rootSplit=gpubytes
branchdp.rootSplit.hostBudgetBytes=280000000000
branchdp.rootSplit.gpuBudgetBytes=23085449216
packstar.admission.gpuWorkPerSecondPerGpu=614000000
packstar.admission.gpuOutOfCoreWorkPerSecondPerGpu=170000000
branchdp.cutoff.strategy=RESIDUAL_BUDGET
packstar.pac.frequencySeverity.jointMomentLearning=true
packstar.pac.frequencySeverity.relativeBoundKcal=1
packstar.pac.frequencySeverity.severityCap=20
packstar.pac.frequencySeverity.severityPremiseId=conditional-relative-gauge-S0-20-not-externally-recalibrated
packstar.pac.frequencySeverity.shrinkGrid=0:0,2:5,5:10,10:20
packstar.pac.frequencySeverity.alphaGrid=1
packstar.pac.frequencySeverity.folds=2
packstar.pac.frequencySeverity.minShiftEssFraction=0.000000000001
packstar.pac.frequencySeverity.discoveryMinShiftEssFraction=0.000000000001
packstar.pac.frequencySeverity.maxRefits=8
packstar.pac.frequencySeverity.discoverySamples=100
packstar.pac.frequencySeverity.discoveryMaxSamples=400
packstar.pac.frequencySeverity.validationSamples=400
packstar.pac.frequencySeverity.tripleEtaSelectionStrategy=budget-forward
packstar.pac.frequencySeverity.tripleEtaMaxWork=1000000000000
packstar.pac.frequencySeverity.tripleEtaMaxMStates=1000000000
packstar.pac.frequencySeverity.tripleEtaMaxTableBytes=64GiB
packstar.pac.frequencySeverity.tripleEtaMaxTotalTableBytes=256GiB
packstar.pac.frequencySeverity.tripleEtaMaxHostBytes=280000000000
packstar.pac.frequencySeverity.tripleEtaMaxFileBytes=512GiB
packstar.pac.frequencySeverity.tripleEtaScale=1
packstar.pac.frequencySeverity.tripleEtaScaleGrid=1
packstar.pac.frequencySeverity.tripleEtaMaxAssignments=1000000
packstar.pac.frequencySeverity.tripleEtaMaxPositionTriples=3
packstar.pac.frequencySeverity.tripleEtaMaxFillEdges=3
packstar.pac.frequencySeverity.tripleEtaMinCellContexts=1
packstar.pac.frequencySeverity.tripleEtaPriorStrength=4
packstar.pac.frequencySeverity.tripleEtaLocalCapKcal=2
packstar.pac.frequencySeverity.tripleEtaResidualCapKcal=3
packstar.pac.frequencySeverity.minTrainCount=5
packstar.pac.frequencySeverity.maxUndertrainedAmplification=1.25
packstar.pac.frequencySeverity.sizeSafety=0.9
packstar.pac.frequencySeverity.severityTestAlpha=0.05
packstar.pac.samples=1000
packstar.pac.trainSamples=500
packstar.pac.pilotSamples=100
packstar.pac.maxEstSamples=4000
packstar.pac.monitorSamples=100
packstar.pac.unreachableCap=400
packstar.pac.nstarInflate=1.3
packstar.pac.confidence=0.05
packstar.pac.targetEpsilon=0.683
packstar.pac.sampling.gpu=false
packstar.pac.sampling.parallel=true
packstar.pac.sampling.threads=48
osprey.wmb.numGpus=0
osprey.wmb.streamsPerGpu=64
osprey.packstar.reduceMinimizations=true
osprey.packstar.correctionTightening=true
osprey.bench.method=packstar_pfunc
osprey.packstarPfunc.maxMut=1
packstar.pac.sampling.timing=true
packstar.pac.dp.timing=true
osprey.bench.epsilon=0.683
osprey.bench.numCPUs=48
'''
    props = dict(line.split('=', 1) for line in fixed.splitlines() if line)
    props.update({'java.io.tmpdir': str(run / 'tmp'),
                  'branchdp.dp.mmap.dir': str(run / 'dp_mmap'),
                  'branchdp.cutoff.residualBudget': case['residual_budget'],
                  'packstar.pac.randomSeed': str(seed),
                  'packstar.pac.frequencySeverity.outputDir': str(run / 'adaptive_frequency_severity'),
                  'osprey.packstarPfunc.state': case['state'],
                  'osprey.packstarPfunc.seqIndex': str(case['seq_index']),
                  'osprey.packstarPfunc.outputTsv': str(run / 'result.tsv'),
                  'osprey.bench.outputDir': str(run), 'osprey.bench.designId': case['design'],
                  'osprey.bench.pdbPath': str(run / 'input.pdb'),
                  'osprey.bench.mutable': case['mutable'], 'osprey.bench.flexible': case['flexible']})
    return props


def audit_arm(run, case, arm):
    rows = table(run / 'result.tsv')
    assert len(rows) == 1, ('expected one result', rows)
    row = rows[0]
    for key, expected in [('design_id', case['design']), ('state', case['state'].capitalize()),
                          ('sequence', case['sequence']), ('seq_index', str(case['seq_index']))]:
        assert row[key].strip() == expected, (key, row[key], expected)
    states = list((run / 'adaptive_frequency_severity').glob('state-*'))
    assert len(states) == 1, states
    state = states[0]
    protocol = keys(state / 'frequency_severity_protocol.tsv')
    assert protocol['proposalLearningEnabled'] == str(arm != 'no-learning').lower()
    assert protocol['tripleEtaEnabled'] == str(arm == 'budget-forward').lower()
    winner = keys(state / 'eta_selected.tsv')
    if winner and arm != 'budget-forward':
        assert winner['tripleEtaSelectedPositionTriples'] == '0', winner
    checked_samples = 0
    if arm == 'no-learning':
        log = (run / 'run.log').read_text(errors='replace')
        assert 'joint-pair-fit' not in log and 'triple-m2-fit' not in log
        if winner:
            assert winner['candidate'] == 'fixed-qm-no-learning', winner
            assert winner['proposalDpSweeps'] == '0', winner
        for path in state.glob('*_samples.tsv'):
            for sample in table(path):
                if 'etaKcal' in sample:
                    assert float(sample['etaKcal']) == 0
                    assert float(sample['eMinKcal']) == float(sample['eProposalKcal'])
                    checked_samples += 1
    if row['status'] == 'Estimated':
        assert math.isfinite(float(row['epsilon'])) and float(row['epsilon']) <= 0.683
    training = table(state / 'adaptive_eta_round_00_frequency_severity_training_samples.tsv')
    fingerprint = [(r['conf'], r['eTrueKcal'], r['eMinKcal']) for r in training]
    evidence = {'ok': True, 'state': str(state), 'status': row['status'],
                'zero_eta_samples_checked': checked_samples,
                'initial_training_count': len(fingerprint),
                'initial_training_sha256': hashlib.sha256(json.dumps(fingerprint).encode()).hexdigest() if fingerprint else None}
    dump(run / 'technical_audit.json', evidence)
    return evidence


def run_task(build, root, task):
    spec = json.loads((build / 'cases.json').read_text())
    case = spec['cases'][task // 3]
    seed = spec['seeds'][task % 3]
    base = root / case['id'] / f'seed_{seed}'
    base.mkdir(parents=True, exist_ok=False)
    assert int(os.environ['SLURM_CPUS_PER_TASK']) == 48
    gpu = subprocess.check_output(['nvidia-smi', '--query-gpu=name,uuid,memory.total', '--format=csv,noheader'], text=True)
    names = gpu.strip().splitlines()
    # Slurm may expose more physical devices on an exclusive eight-GPU node.
    # CUDA_VISIBLE_DEVICES constrains this run; DP maxGpus is explicitly four.
    visible = os.environ.get('CUDA_VISIBLE_DEVICES', '').split(',')
    assert len(visible) >= 4 and all('A5000' in line for line in names), (visible, gpu)
    allocated_visible = ','.join(visible)
    # Exclusive gres allocations on fennario expose all eight GPUs even when
    # four were requested. Restrict this JVM to four allocated A5000 devices.
    os.environ['CUDA_VISIBLE_DEVICES'] = ','.join(visible[:4])
    (base / 'gpus.csv').write_text(gpu)
    checked_run(['lscpu'], base / 'lscpu.txt')
    checked_run(['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID']], base / 'slurm_job.txt')
    order = spec['arms'][task % 3:] + spec['arms'][:task % 3]
    metadata = {'case': case, 'seed': seed, 'arm_order': order, 'build': str(build),
                'git_head': (build / 'git_head.txt').read_text().strip(),
                'node': os.environ.get('SLURMD_NODENAME'), 'partition': os.environ.get('SLURM_JOB_PARTITION'),
                'job': os.environ['SLURM_JOB_ID'], 'task': task, 'gpu_inventory': gpu,
                'cuda_visible_devices': os.environ.get('CUDA_VISIBLE_DEVICES'),
                'allocated_cuda_visible_devices': allocated_visible,
                'account': os.environ.get('SLURM_JOB_ACCOUNT')}
    dump(base / 'task.json', metadata)
    input_manifest = json.loads((build / 'input_manifest.json').read_text())
    technical_ok = True
    audits = {}
    for arm_index, arm in enumerate(order):
        run = base / arm
        cache = run / 'emat_cache' / case['design']
        cache.mkdir(parents=True)
        for folder in ('tmp', 'cuda_cache', 'dp_mmap'):
            (run / folder).mkdir()
        copied = {}
        for name in ['input.pdb'] + [p.name for p in (build / 'inputs' / case['design']).glob('*.dat')]:
            source = build / 'inputs' / case['design'] / name
            target = run / name if name == 'input.pdb' else cache / name
            subprocess.run(['cp', '--reflink=auto', str(source), str(target)], check=True)
            digest = sha(target)
            assert digest == input_manifest[case['design'] + '/' + name]['sha256']
            copied[name] = digest
        dump(run / 'input_hashes.json', copied)
        props = common_properties(case, seed, run)
        props['packstar.pac.frequencySeverity.proposalLearning'] = str(arm != 'no-learning').lower()
        props['packstar.pac.frequencySeverity.tripleEta'] = str(arm == 'budget-forward').lower()
        dump(run / 'properties.json', props)
        command = java_command(build, props)
        (run / 'command.sh').write_text(shlex.join(command) + '\n')
        env = dict(os.environ, TMPDIR=str(run / 'tmp'), CUDA_CACHE_PATH=str(run / 'cuda_cache'))
        meta = dict(metadata, arm=arm, arm_order_index=arm_index)
        dump(run / 'run.json', meta)
        print(f'START case={case["id"]} seed={seed} arm={arm} node={metadata["node"]}', flush=True)
        with (run / 'run.log').open('w') as stream:
            proc = subprocess.run(['/usr/bin/time', '-o', str(run / 'wall.time'),
                                   '-f', 'elapsed=%e maxRssKiB=%M', 'timeout',
                                   '--signal=TERM', '--kill-after=60s', '5h', *command],
                                  cwd=build / 'source', env=env,
                                  stdout=stream, stderr=subprocess.STDOUT)
        meta['java_exit'] = proc.returncode
        dump(run / 'run.json', meta)
        if proc.returncode == 0:
            try:
                audits[arm] = audit_arm(run, case, arm)
            except Exception as exc:
                technical_ok = False
                dump(run / 'technical_audit.json', {'ok': False, 'error': repr(exc)})
        elif proc.returncode not in (124, 137):
            technical_ok = False
        print(f'END case={case["id"]} seed={seed} arm={arm} exit={proc.returncode}', flush=True)
        # Only this invocation's generated workspace; keep all scientific audits.
        if proc.returncode == 0:
            shutil.rmtree(run / 'dp_mmap')
    fingerprints = {a['initial_training_sha256'] for a in audits.values() if a['initial_training_sha256']}
    counts = {a['initial_training_count'] for a in audits.values() if a['initial_training_sha256']}
    if len(fingerprints) > 1 or (counts and counts != {500}):
        technical_ok = False
    dump(base / 'paired_audit.json', {'ok': technical_ok, 'arms_audited': list(audits),
         'initial_training_fingerprints_equal': len(fingerprints) == 1,
         'initial_training_counts': sorted(counts), 'audits': audits})
    # Gate requires all three complete technical audits. Statistical failure is valid.
    if task == 0 and len(audits) != 3:
        technical_ok = False
    if not technical_ok:
        raise SystemExit('Technical audit failed; inspect preserved outputs.')


def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def summarize(build, root):
    root.mkdir(parents=True, exist_ok=True)
    spec = json.loads((build / 'cases.json').read_text())
    rows = []
    for case in spec['cases']:
        for seed in spec['seeds']:
            for arm in spec['arms']:
                run = root / case['id'] / f'seed_{seed}' / arm
                meta = json.loads((run / 'run.json').read_text()) if (run / 'run.json').exists() else {}
                results = table(run / 'result.tsv')
                result = results[0] if results else {}
                states = list((run / 'adaptive_frequency_severity').glob('state-*'))
                state = states[0] if len(states) == 1 else run / '_missing_state'
                final = keys(state / 'frequency_severity_final_interval.tsv')
                winner = keys(state / 'eta_selected.tsv')
                failure = keys(state / 'frequency_severity_failure.tsv')
                audit = json.loads((run / 'technical_audit.json').read_text()) if (run / 'technical_audit.json').exists() else {}
                status = result.get('status', 'MISSING')
                exit_code = meta.get('java_exit')
                if exit_code in (124, 137):
                    status = 'TIMEOUT'
                elif exit_code not in (None, 0):
                    status = 'PROCESS_ERROR'
                elif audit.get('ok') is False:
                    status = 'AUDIT_FAILED'
                elif run.exists() and not results:
                    status = 'INCOMPLETE'
                wall_text = (run / 'wall.time').read_text() if (run / 'wall.time').exists() else ''
                wall = re.search(r'elapsed=([\d.]+) maxRssKiB=(\d+)', wall_text)
                log = (run / 'run.log').read_text(errors='replace') if (run / 'run.log').exists() else ''
                totals = re.findall(r'\[PACK\*\] Total: (\d+) ms, (\d+) unique CCD calls for (\d+) sample records', log)
                rounds = [keys(p) for p in sorted(state.glob('adaptive_eta_round_*_selection.tsv'))]
                rows.append(dict(case=case['id'], design=case['design'], state=case['state'],
                    sequence=case['sequence'], seq_index=case['seq_index'], rb=case['residual_budget'],
                    seed=seed, arm=arm, status=status, raw_status=result.get('status'),
                    technical_ok=audit.get('ok'), node=meta.get('node'), partition=meta.get('partition'),
                    wall_s=float(wall[1]) if wall else None, peak_rss_kib=int(wall[2]) if wall else None,
                    pfunc_s=numeric(result.get('elapsed_s')), epsilon=numeric(result.get('epsilon')),
                    lower_log10=numeric(result.get('lower_log10')), upper_log10=numeric(result.get('upper_log10')),
                    selected_k=numeric(final.get('selectedTripleEtaPositionTriples', winner.get('tripleEtaSelectedPositionTriples'))),
                    selected_scopes=final.get('selectedTripleEtaPositionScopes', winner.get('tripleEtaSelectedPositionScopes')),
                    selected_candidate=final.get('selectedCandidate', winner.get('candidate')),
                    final_n=numeric(final.get('finalCcd', failure.get('finalCcd'))),
                    unique_ccd=numeric(final.get('totalCcd', failure.get('totalCcd'))),
                    logical_samples=int(totals[-1][2]) if totals else numeric(failure.get('totalSampleRecords')),
                    fitting_rounds=len(rounds), rounds_with_selected_triples=sum(float(r.get('selectedTripleEtaPositionTriples', 0)) > 0 for r in rounds),
                    corrected_dp_sweeps=numeric(winner.get('proposalDpSweeps')),
                    failure_reason=failure.get('reason'), java_exit=exit_code, run=str(run)))
    with (root / 'comparison.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    paired = []
    for case in spec['cases']:
        for seed in spec['seeds']:
            group = {r['arm']: r for r in rows if r['case'] == case['id'] and r['seed'] == seed}
            pair, triple = group['pair-only'], group['budget-forward']
            paired.append({'case': case['id'], 'seed': seed, 'pair_status': pair['status'],
                           'triple_status': triple['status'], 'triple_selected_k': triple['selected_k'],
                           'triple_only_success': triple['status'] == 'Estimated' and pair['status'] == 'Aborted',
                           'pair_only_success': pair['status'] == 'Estimated' and triple['status'] == 'Aborted',
                           'wall_triple_over_pair': triple['wall_s'] / pair['wall_s'] if pair['status'] == triple['status'] == 'Estimated' and pair['wall_s'] and triple['wall_s'] else None})
    counts = {arm: {s: sum(r['arm'] == arm and r['status'] == s for r in rows)
                    for s in sorted({r['status'] for r in rows if r['arm'] == arm})} for arm in spec['arms']}
    dump(root / 'comparison.json', {'planned_runs': len(rows), 'status_counts': counts, 'paired': paired, 'rows': rows})
    lines = ['# Current-version hard-case ablation', '',
             f'{len(spec["cases"])} predeclared sequence-state case(s) × three seeds × three arms. RB=1 throughout. Same frozen build and target within each case; 48 CPUs/four A5000, same node within each seed. Total wall includes model selection and all CCD. Historical failed-then-successful-with-final-triples observations selected the cases; they are not matched controls.', '',
             f'| Arm | Estimated / {len(spec["cases"]) * len(spec["seeds"])} | All outcomes |', '|---|---:|---|']
    for arm in spec['arms']:
        lines.append(f'| {arm} | {counts[arm].get("Estimated", 0)} | {counts[arm]} |')
    lines += ['', '| Case | Seed | Pair | Optional triple | K | Triple/pair wall on joint successes |', '|---|---:|---|---|---:|---:|']
    for pair in paired:
        ratio = pair['wall_triple_over_pair']
        lines.append('| ' + ' | '.join(str(x) for x in [pair['case'], pair['seed'], pair['pair_status'], pair['triple_status'], pair['triple_selected_k'], f'{ratio:.3f}' if ratio is not None else 'NA']) + ' |')
    lines += ['', 'Triple-only success means a valid Estimated triple-arm result paired with a statistical Aborted pair-only result. Timeouts, audit errors and missing runs are separately counted. K=0 is allowed. Selected K>0 alone does not prove benefit. S0=20 remains a conditional premise; this experiment does not recalibrate it.', '',
              'Full per-run time, CCD, final N, intervals, failure reasons and output paths: comparison.tsv. Initial training is checked across arms; later adaptive data may differ.']
    (root / 'comparison.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(counts, indent=2))


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Submit through Slurm with account=grisman.')
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['prepare', 'run', 'summarize'])
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--task', type=int)
    args = parser.parse_args()
    storage = Path('/usr/xtmp/lz280').resolve()
    assert storage in args.build.resolve().parents
    if args.root:
        assert storage in args.root.resolve().parents
    if args.mode == 'prepare':
        prepare(args.build)
    elif args.mode == 'run':
        run_task(args.build, args.root, args.task)
    else:
        summarize(args.build, args.root)


if __name__ == '__main__':
    main()
