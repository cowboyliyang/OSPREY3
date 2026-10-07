"""Submit exactly five approved additions without increasing local lane counts."""
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from collections import Counter
from datetime import datetime
from zoneinfo import ZoneInfo

assert os.environ.get('SLURM_JOB_ID'), 'Submit through Slurm with account grisman.'
HERE = Path(__file__).resolve().parent
FOCUS = HERE.parent/'deadline_focus_20261006'
INVENTORY = HERE.parent/'current_inventory_20261006'
PLAN = json.loads((HERE/'plan.json').read_text())
record = json.loads((HERE/'registry.json').read_text())
assert record['status']=='PREPARED' and not record['jobs'], 'Inspect existing submissions instead of repeating.'
prep = Path(record['prep'])
assert (prep/'READY').is_file()
audit = prep/'submission'
audit.mkdir(exist_ok=False)
commands = []

def save(path, value):
    temporary = path.with_name(path.name+'.next')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')
    temporary.replace(path)

def checkpoint():
    save(HERE/'registry.json', record)
    save(audit/'registry.json', record)

def run(*args):
    result = subprocess.run(args, capture_output=True, text=True, timeout=60)
    commands.append(dict(arguments=list(args), exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr))
    save(audit/'commands.json', commands)
    result.check_returncode()
    return result.stdout

def state(job):
    return dict(re.findall(r'(?:^|\s)([A-Za-z][A-Za-z0-9_]*)=(\S+)', run('scontrol', 'show', 'job', job, '-o')))

def never_started(fields):
    assert fields['JobState']=='PENDING' and fields['RunTime']=='00:00:00' and fields.get('Restarts')=='0', fields
    assert fields['Account']=='grisman' and fields['NumCPUs']=='64', fields

def parent(text):
    if text in ('', '(null)', '0'):
        return None
    match = re.fullmatch(r'afterany:(\d+(?:_\d+)?)(?:\((?:un)?fulfilled\))?', text)
    assert match, text
    return match[1]

def queue():
    raw = run('squeue', '-h', '-r', '-u', 'lz280', '-o', '%i|%j|%T|%P|%E|%C|%a')
    rows = {}
    for line in raw.splitlines():
        parts = line.split('|')
        if len(parts)==7 and parts[1].startswith('markstar') and parts[5]=='64':
            rows[parts[0]] = dict(zip(('name', 'state', 'partition', 'dependency', 'cpus', 'account'), parts[1:]))
    return raw, rows

def validate(rows, links):
    assert all(row['account']=='grisman' for row in rows.values())
    children = Counter(p for p in links.values() if p in rows)
    assert all(n<=1 for n in children.values()), ('A serial lane would branch', children)
    roots = Counter(rows[j]['partition'] for j, p in links.items() if p not in rows)
    assert roots['compsci']<=PLAN['local_limits']['compsci'], roots
    assert roots['grisman']<=PLAN['local_limits']['fennario'], roots
    for job, predecessor in links.items():
        if predecessor in rows:
            assert rows[job]['partition']==rows[predecessor]['partition'], (job, predecessor)
        seen, cursor = set(), job
        while cursor in rows:
            assert cursor not in seen, 'Dependency cycle'
            seen.add(cursor)
            cursor = links[cursor]
    return dict(roots)

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    temporary = path.with_name(path.name+'.next')
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)

raw, before = queue()
(audit/'queue_before.txt').write_text(raw)
links = {job: parent(row['dependency']) for job, row in before.items()}
record['lane_roots_before'] = validate(before, links)
table = {r['design']: r for r in read(INVENTORY/'current_plan.tsv')}
definitions = {r['design_id']: r for r in read(prep/'package/designs.tsv')}
requests = [dict(r) for r in record['requests']]
child_change = None
for request in requests:
    name = request['design']
    assert table[name]['state']=='待提交' and table[name]['recommendation']=='保留' and not table[name]['job']
    expected_partition = 'compsci' if request['group']=='compsci' else 'grisman'
    predecessor = state(request['after_job'])
    assert predecessor['Account']=='grisman' and predecessor['NumCPUs']=='64'
    assert predecessor['Partition']==expected_partition
    if request.get('before_job'):
        fields = state(request['before_job'])
        if fields['JobState']=='PENDING':
            never_started(fields)
            assert parent(fields['Dependency'])==request['after_job'] and fields['Reason']!='JobHeldUser'
            child_change = dict(job=request['before_job'], old_parent=request['after_job'], design=name, before=fields)
        else:
            # Preserve any work which began after the planning snapshot.
            assert fields['JobState'] in ('RUNNING', 'COMPLETING', 'COMPLETED'), fields
            request['after_job'] = request.pop('before_job')
            request['after_design'] = request.pop('before_design')
            request['adaptation'] = 'Append after p6 because it has already started.'
    children = [job for job, predecessor in links.items() if predecessor==request['after_job']]
    assert children==([request['before_job']] if request.get('before_job') else []), (request, children)

planned_rows = {j: dict(r) for j, r in before.items()}
planned_links = dict(links)
for request in requests:
    placeholder = 'new:'+request['design']
    planned_rows[placeholder] = dict(account='grisman', partition='compsci' if request['group']=='compsci' else 'grisman')
    planned_links[placeholder] = request['after_job']
    if request.get('before_job'):
        planned_links[request['before_job']] = placeholder
record['planned_lane_roots'] = validate(planned_rows, planned_links)
protected = {job: state(job) for job, row in before.items() if row['state']=='RUNNING'}
save(audit/'protected_running_before.json', protected)
for path in (FOCUS/'design_actions.tsv', FOCUS/'local_lane_dependencies.tsv', FOCUS/'policy.json',
             FOCUS/'registry.json', INVENTORY/'reporting_preferences.json'):
    backup = audit/'before'/path.parent.name/path.name
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup)
record.update(status='SUBMITTING', effective_requests=requests, audit=str(audit), child_change=child_change)
checkpoint()
try:
    if child_change:
        never_started(state(child_change['job']))
        run('scontrol', 'hold', child_change['job'])
        never_started(state(child_change['job']))
        child_change['held_by_management'] = True
        checkpoint()
    for request in requests:
        name, group = request['design'], request['group']
        index = str(definitions[name]['group_index'])
        options = ['--partition=compsci', '--exclude=linux[31-40]'] if group=='compsci' else [
            '--partition=grisman', '--constraint=a5000', '--exclude=grisman-37,grisman-40,jerry[1-7]']
        receipt = run('sbatch', '--parsable', '--account=grisman', '--hold',
            '--job-name=markstar-add-'+name.replace('_flex_', '-'), *options,
            '--array='+index+'-'+index+'%1', '--dependency=afterany:'+request['after_job'],
            '--export=ALL,MARKSTAR_EXTENSION_PREP='+str(prep)+',MARKSTAR_EXTENSION_GROUP='+group,
            str(prep/'markstar_array.slurm')).strip().split(';')[0]
        assert receipt.isdigit(), receipt
        job = receipt+'_'+index
        entry = dict(design=name, job=job, array_job=receipt, group=group,
                     after_job=request['after_job'], after_design=request['after_design'],
                     prep=str(prep), output=str(prep.parent/('markstar_A'+receipt)/name),
                     cpus=64, memory_gib=192, limit_hours=336)
        record['jobs'].append(entry)
        checkpoint()
        info = state(job)
        never_started(info)
        assert parent(info['Dependency'])==request['after_job']
        assert info['MinMemoryNode']=='192G' and info['TimeLimit']=='14-00:00:00', info
    by_name = {r['design']: r for r in record['jobs']}
    if child_change:
        never_started(state(child_change['job']))
        child_change['new_parent'] = by_name[child_change['design']]['job']
        run('scontrol', 'update', 'JobId='+child_change['job'], 'Dependency=afterany:'+child_change['new_parent'])
        assert parent(state(child_change['job'])['Dependency'])==child_change['new_parent']
        checkpoint()
    raw, after = queue()
    (audit/'queue_submitted_held.txt').write_text(raw)
    record['lane_roots_after_submission'] = validate(after, {j: parent(r['dependency']) for j, r in after.items()})
    edges = read(FOCUS/'local_lane_dependencies.tsv')
    if child_change:
        edge = next(r for r in edges if r['job']==child_change['job'])
        edge.update(parent=child_change['new_parent'], action='COVERAGE_SMALLER_TIER_FIRST_20261007')
    edges.extend(dict(job=r['job'], parent=r['after_job'], action='QUEUED_APPROVED_PLAN_ADDITION') for r in record['jobs'])
    write(FOCUS/'local_lane_dependencies.tsv', edges)
    actions = read(FOCUS/'design_actions.tsv')
    for entry in record['jobs']:
        row = next(r for r in actions if r['design']==entry['design'])
        row.update(action='KEEP_USER_APPROVED_ADDITION', current_assignment='local', proposed_destination=entry['group'],
                   current_local_job=entry['job'], current_status='SUBMITTED_CPU_ONLY')
    write(FOCUS/'design_actions.tsv', actions)
    write(HERE/'local_submissions.tsv', record['jobs'])
    prefs = json.loads((INVENTORY/'reporting_preferences.json').read_text())
    for addition in prefs['plan_additions']:
        if addition['design'] in by_name:
            addition['destination'] = by_name[addition['design']]['group']
            addition['submission_registry'] = str(HERE/'registry.json')
    prefs['submitted_additions_registry'] = str(HERE/'registry.json')
    prefs['scope'] = 'Approved unified plan; all five local additions are submitted. Conditional tiers remain table-only; DCC receipt gaps remain explicit.'
    save(INVENTORY/'reporting_preferences.json', prefs)
    for path in (FOCUS/'policy.json', FOCUS/'registry.json'):
        content = json.loads(path.read_text())
        content['queued_five_registry'] = str(HERE/'registry.json')
        save(path, content)
    record['status'] = 'SUBMITTED_HELD'
    checkpoint()
    for entry in record['jobs']:
        run('scontrol', 'release', entry['job'])
    if child_change:
        run('scontrol', 'release', child_change['job'])
    raw, final_queue = queue()
    (audit/'queue_after.txt').write_text(raw)
    record['lane_roots_after'] = validate(final_queue, {j: parent(r['dependency']) for j, r in final_queue.items()})
    for entry in record['jobs']:
        info = state(entry['job'])
        assert info['Reason']!='JobHeldUser' and info['Account']=='grisman' and info['NumCPUs']=='64', info
        if info['JobState']=='PENDING':
            assert parent(info['Dependency']) in (entry['after_job'], None), info
        entry['state_after_submission'] = info['JobState']
    if child_change:
        assert state(child_change['job'])['Reason']!='JobHeldUser'
    for job, old in protected.items():
        current = state(job)
        assert all(current[key]==old[key] for key in ('StartTime', 'Command', 'NumCPUs', 'NodeList', 'TimeLimit')), job
    record.update(status='SUBMITTED', applied=datetime.now(ZoneInfo('America/New_York')).isoformat(),
                  protected_running_tasks=len(protected), conditional_jobs_restarted=[], dcc_changes=[], new_packstar_jobs=[])
    checkpoint()
except BaseException as error:
    record.update(status='NEEDS_REVIEW', failure=repr(error))
    checkpoint()
    raise
print('SUBMITTED_FIVE', json.dumps(record, ensure_ascii=False), flush=True)
