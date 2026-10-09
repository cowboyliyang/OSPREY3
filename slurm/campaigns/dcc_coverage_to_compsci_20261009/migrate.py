"""Preflight, submit and register four approved DCC-to-compsci coverage runs."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from collections import Counter
from datetime import datetime
from zoneinfo import ZoneInfo

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm under account grisman.'
HERE = Path(__file__).resolve().parent
INV = HERE.parent / 'current_inventory_20261006'
FOCUS = HERE.parent / 'deadline_focus_20261006'
NAMES = ['3bua_flex_m4', '3ma2_flex_p0', '4wwi_flex_p1', '4wyu_flex_p0']
ROOT = Path('/usr/xtmp/lz280/markstar_dcc_coverage_to_compsci_20261009')
OUT = ROOT / ('prep_' + os.environ['SLURM_JOB_ID'])
assert not (HERE / 'registry.json').exists(), 'Inspect existing receipts before repeating migration.'
BASE = Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
NEW = Path('/usr/xtmp/lz280/markstar_deadline_focus_20261006/prep_12826915/package')

def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def save(path, value):
    temporary = path.with_name(path.name + '.next')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)

def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    temporary = path.with_name(path.name + '.next')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

definitions = {}
for source in (BASE, NEW):
    for row in read(source / 'designs.tsv'):
        definitions[row['design_id']] = (row, source)
historical = {r['design']: r for r in read(Path('/usr/xtmp/lz280/markstar_extension_inventory_12826374/all161.tsv'))}
table = {r['design']: r for r in read(INV / 'current_plan.tsv')}
for name in NAMES:
    assert table[name]['destination'] == 'dcc' and table[name]['state'].startswith('DCC'), table[name]
    assert name not in json.loads((INV / 'local_job_overrides.json').read_text())
sources = [source / row['pdb_relative'] for name in NAMES for row, source in [definitions[name]]]
print('STORAGE_ESTIMATE', json.dumps(dict(pdb_files=4, pdb_bytes=sum(p.stat().st_size for p in sources),
      preparation_files_under=80, preparation_bytes_under=50 * 2**20,
      runtime='Four MARK runs, existing compiled classpaths reused; outputs and caches stay under scratch.')), flush=True)
OUT.mkdir(parents=True, exist_ok=False)
package = OUT / 'package'
(package / 'structures').mkdir(parents=True)
(package / 'preflight').mkdir()
shutil.copy2(__file__, OUT / 'migrate.py')
commands = []

def run(*args):
    result = subprocess.run(args, capture_output=True, text=True, timeout=60)
    commands.append(dict(arguments=list(args), exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr))
    save(OUT / 'scheduler_commands.json', commands)
    result.check_returncode()
    return result.stdout.strip()

def info(job):
    return dict(re.findall(r'(?:^|\s)(\w+)=(\S+)', run('scontrol', 'show', 'job', job, '-o')))

packs = {}
for group in ('packstar_gpu101_20261005/run_12814879', 'markstar_deadline_focus_20261006/pack_12827802'):
    for row in read(Path('/usr/xtmp/lz280') / group / 'status.tsv'):
        if row['design'] in NAMES and row['status'] in ('COMPLETED', 'INCOMPLETE_ESTIMATES'):
            packs[row['design']] = row
assert set(packs) == set(NAMES), packs
rows, configs, verification = [], {}, []
for index, name in enumerate(NAMES):
    original, source_package = definitions[name]
    row = dict(original)
    source_pdb = source_package / row['pdb_relative']
    assert digest(source_pdb) == row['pdb_sha256']
    pdb = package / row['pdb_relative']
    shutil.copy2(source_pdb, pdb)
    assert digest(pdb) == row['pdb_sha256']
    anchor = name if name in historical else row['anchor_design']
    source = Path(historical[anchor]['mark_path'])
    command = shlex.split((source / 'command.sh').read_text())
    updates = {'-XX:ActiveProcessorCount=': '64', '-Dosprey.bench.numCPUs=': '64',
               '-Dosprey.bench.designId=': name, '-Dosprey.bench.outputDir=': '{RUN}',
               '-Djava.io.tmpdir=': '{RUN}/tmp', '-Dosprey.bench.pdbPath=': str(pdb),
               '-Dosprey.bench.mutable=': row['mutable'], '-Dosprey.bench.flexible=': row['flexible']}
    rewritten, changed = [], set()
    for token in command:
        key = next((k for k in updates if token.startswith(k)), None)
        if key:
            token = key + updates[key]
            changed.add(key)
        rewritten.append(token.replace(str(source), '{RUN}'))
    assert changed == set(updates) and '-Dosprey.bench.method=markstar' in rewritten
    assert all(Path(p).exists() for p in rewritten[rewritten.index('-cp') + 1].split(':'))
    heap = int(re.fullmatch(r'-Xmx(\d+)g', next(t for t in rewritten if t.startswith('-Xmx')))[1])
    assert 0 < heap < 192
    pack_root = Path(packs[name]['output'])
    pack_manifest = json.loads((pack_root / 'manifest.json').read_text())
    for key in ('mutable', 'flexible', 'pdb_sha256', 'expected_sequences'):
        assert str(pack_manifest['design'][key]) == str(row[key]), (name, key)
    pack_results = list(csv.DictReader((pack_root / (name + '_packstar.csv')).open()))
    assert len(pack_results) == int(row['expected_sequences'])
    preflight = package / 'preflight' / name
    preflight.mkdir()
    (preflight / 'tmp').mkdir()
    precommand = []
    for token in rewritten:
        token = token.replace('{RUN}', str(preflight))
        if token.startswith('-Xmx'):
            token = '-Xmx8g'
        elif token.startswith('-XX:ActiveProcessorCount='):
            token = '-XX:ActiveProcessorCount=4'
        elif token.startswith('-Dosprey.bench.numCPUs='):
            token = '-Dosprey.bench.numCPUs=4'
        elif token.startswith('-Dosprey.bench.method='):
            token = '-Dosprey.bench.method=sequence_dump'
        precommand.append(token)
    (preflight / 'command.sh').write_text(shlex.join(precommand) + '\n')
    with (preflight / 'run.log').open('w') as log:
        subprocess.run(precommand, cwd=preflight, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    log = (preflight / 'run.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log
    positions = re.findall(r'\[FRONTIER_POSITION\] state=Complex residue=(\S+)', log)
    assert len(positions) == int(row['total_positions'])
    assert set(positions) == set((row['mutable'] + ';' + row['flexible']).split(';'))
    sequences = read(preflight / (name + '_sequences.tsv'))
    assert len(sequences) == int(row['expected_sequences'])
    assert {r['sequence'] for r in sequences} == {r['sequence'] for r in pack_results}
    row.update(destination='compsci', group_index=index, requested_cpus=64, requested_memory_gib=192,
               java_heap_gib=heap, limit_hours=336, kind='USER_APPROVED_DCC_COVERAGE_MIGRATION',
               reference_design=anchor, reference_mark_path=str(source), source_package=str(source_package))
    rows.append(row)
    configs[name] = dict(row=row, command_template=rewritten, input_pdb=str(pdb), model_contains='AMD EPYC 9554',
                         source_command_sha256=digest(source / 'command.sh'))
    check = dict(design=name, input_preflight='PASS', sequence_set_matches_PACK=True,
                 expected_sequences=len(sequences), positions=len(positions), heap_gib=heap,
                 historical_completed_under14days=name in historical and historical[name]['kind'] == 'completed'
                 and float(historical[name]['mark_hours']) < 336)
    verification.append(check)
    print('PREFLIGHT_OK', json.dumps(check), flush=True)
write(package / 'designs.tsv', rows)
save(package / 'run_configs.json', configs)
for filename in ('markstar_runner.py', 'markstar_array.slurm'):
    if filename.endswith('.py'):
        compile((HERE / filename).read_text(), filename, 'exec')
    else:
        subprocess.run(['bash', '-n', str(HERE / filename)], check=True)
    shutil.copy2(HERE / filename, OUT / filename)
(OUT / 'READY').write_text('Four input preflights passed; sequence sets match measured PACK results.\n')
write(HERE / 'designs.tsv', rows)
record = dict(status='PREPARED', management_job=os.environ['SLURM_JOB_ID'], prep=str(OUT), package=str(package),
              verification=verification, jobs=[], remote_scheduler_accessed=False, remote_state='UNKNOWN',
              cpu_model='AMD EPYC 9554', hardware_note='Separate compsci cohort; do not assume DCC CPU matches.')

def checkpoint():
    save(HERE / 'registry.json', record)
    save(OUT / 'registry.json', record)

checkpoint()
receipt = run('sbatch', '--parsable', '--hold', '--account=grisman', '--partition=compsci',
              '--job-name=markstar-dcc-coverage4', '--array=0-3%4', '--exclude=linux[31-40]',
              '--export=ALL,MARKSTAR_EXTENSION_PREP=' + str(OUT), str(OUT / 'markstar_array.slurm')).split(';')[0]
assert receipt.isdigit(), receipt
for index, name in enumerate(NAMES):
    job = receipt + '_' + str(index)
    record['jobs'].append(dict(design=name, job=job, array_job=receipt, destination='compsci', cpus=64,
          memory_gib=192, time_limit_hours=336, previous_destination='dcc', remote_state='UNKNOWN',
          prep=str(OUT), output=str(ROOT / ('markstar_A' + receipt) / name)))
checkpoint()
for entry in record['jobs']:
    state = info(entry['job'])
    assert state['JobState'] == 'PENDING' and state['Reason'] == 'JobHeldUser'
    assert state['Account'] == 'grisman' and state['NumCPUs'] == '64' and state['MinMemoryNode'] == '192G'
    assert state['Partition'] == 'compsci' and state['TimeLimit'] == '14-00:00:00'
    assert state['Dependency'] == '(null)'

# Replace only the admission controller, carrying its existing held backlog forward.
prefs = json.loads((INV / 'reporting_preferences.json').read_text())
old_controller = prefs['queue_backlog_policy']['controller_job']
old_state = info(old_controller)
assert old_state['JobName'] == 'mark-backfill-controller' and old_state['NumCPUs'] == '1'
assert old_state['Account'] == 'grisman' and old_state['JobState'] == 'RUNNING'
old_status = json.loads(Path(prefs['queue_backlog_policy']['controller_status']).read_text())
assert old_status['controller'] == old_controller
managed = [dict(job=j, pool='grisman' if j in ('12853420_0', '12853421_3') else 'compsci',
                adopt_hold=j in old_status['waiting']) for j in prefs['queue_backlog_policy']['managed_jobs']]
managed.extend(dict(job=e['job'], pool='compsci', adopt_hold=True) for e in record['jobs'])
config = dict(limits={'compsci': 10, 'grisman': 5}, jobs=managed, poll_seconds=30,
              replaced_controller=old_controller, migration_registry=str(HERE / 'registry.json'))
save(INV / 'free_backfill_jobs.json', config)
for filename in ('free_backfill.py', 'free_backfill.slurm'):
    shutil.copy2(INV / filename, OUT / filename)
subprocess.run(['bash', '-n', str(OUT / 'free_backfill.slurm')], check=True)
compile((OUT / 'free_backfill.py').read_text(), 'free_backfill.py', 'exec')
controller = run('sbatch', '--parsable', '--account=grisman', '--dependency=afterany:' + old_controller,
                 '--export=ALL,MARKSTAR_BACKFILL_CONFIG=' + str(INV / 'free_backfill_jobs.json') +
                 ',MARKSTAR_BACKFILL_SCRIPT=' + str(OUT / 'free_backfill.py'),
                 str(OUT / 'free_backfill.slurm')).split(';')[0]
assert controller.isdigit()
record.update(status='SUBMITTED_HELD', replaced_controller=old_controller, controller_job=controller)
checkpoint()

for path in (INV / 'local_job_overrides.json', INV / 'reporting_preferences.json', FOCUS / 'policy.json',
             FOCUS / 'design_actions.tsv', FOCUS / 'dcc_actions.tsv', FOCUS / 'dcc_launch_priority.tsv',
             FOCUS / 'dcc_new_cases.txt', FOCUS / 'dcc_new_designs.tsv', FOCUS / 'dcc_keep_existing_pending_cases.txt'):
    backup = OUT / 'before' / path.parent.name / path.name
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup)
overrides = json.loads((INV / 'local_job_overrides.json').read_text())
overrides.update({e['design']: e for e in record['jobs']})
save(INV / 'local_job_overrides.json', overrides)
prefs['updated_on'] = '2026-10-09'
prefs['dcc_coverage_migration_20261009'] = record['jobs']
prefs['scope'] = '89 retained designs; four DCC coverage designs submitted on compsci; remaining DCC states unverified.'
prefs['queue_backlog_policy'].update(controller_job=controller, managed_jobs=[e['job'] for e in managed],
        configuration=str(INV / 'free_backfill_jobs.json'), requested_on='2026-10-09',
        admission='One controller releases explicitly enrolled held tasks when the pool has capacity; prior backlog retained before the four migrated coverage tasks.')
save(INV / 'reporting_preferences.json', prefs)
policy = json.loads((FOCUS / 'policy.json').read_text())
policy['slots']['compsci'] = 10
for section in ('coverage_pending', 'coverage_new'):
    for name in NAMES:
        if name in policy[section]:
            policy[section][name] = 'compsci'
policy['dcc_coverage_migration_20261009'] = dict(registry=str(HERE / 'registry.json'), designs=NAMES,
    remote_state='UNKNOWN', remote_scheduler_accessed=False, remove_from_DCC_launch_lists=True)
save(FOCUS / 'policy.json', policy)
actions = read(FOCUS / 'design_actions.tsv')
by_name = {e['design']: e for e in record['jobs']}
for row in actions:
    if row['design'] in by_name:
        row.update(action='KEEP_USER_MIGRATED_COVERAGE', current_assignment='local', proposed_destination='compsci',
                   current_local_job=by_name[row['design']]['job'], current_status='SUBMITTED_COMPSCI_REMOTE_UNVERIFIED')
write(FOCUS / 'design_actions.tsv', actions)
dcc_actions = read(FOCUS / 'dcc_actions.tsv')
for row in dcc_actions:
    if row['design_id'] in NAMES:
        row.update(action='MOVED_TO_COMPSCI_DO_NOT_LAUNCH', reported_status='REMOTE_STATE_UNVERIFIED_LOCAL_SUBMITTED',
                   cpu_rule='User moved this tier to AMD EPYC 9554 compsci; preserve separate hardware provenance.')
write(FOCUS / 'dcc_actions.tsv', dcc_actions)
policy['dcc_counts'] = dict(Counter(r['action'] for r in dcc_actions))
save(FOCUS / 'policy.json', policy)
for filename in ('dcc_launch_priority.tsv', 'dcc_new_designs.tsv'):
    write(FOCUS / filename, [r for r in read(FOCUS / filename) if r['design_id'] not in NAMES])
for filename in ('dcc_new_cases.txt', 'dcc_keep_existing_pending_cases.txt'):
    path = FOCUS / filename
    path.write_text('\n'.join(n for n in path.read_text().splitlines() if n not in NAMES) + '\n')
save(HERE / 'dcc_handoff.json', dict(date='2026-10-09', designs=NAMES, action='MOVED_TO_COMPSCI_DO_NOT_LAUNCH',
     remote_scheduler_accessed=False, remote_state='UNKNOWN', local_jobs=record['jobs'],
     note='No DCC job was queried or cancelled. Reconcile any existing remote copies on DCC using this migration receipt.'))
run('scancel', old_controller)
record.update(status='SUBMITTED', applied=datetime.now(ZoneInfo('America/New_York')).isoformat(),
              dcc_launch_lists_updated=True, admission_limits=config['limits'])
checkpoint()
print('MIGRATED_FOUR', json.dumps(dict(jobs=record['jobs'], controller_job=controller), ensure_ascii=False), flush=True)
