"""Remove only the three omitted pending MARK* jobs and preserve CPU lanes."""
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime
from zoneinfo import ZoneInfo

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm under account grisman.'
HERE = Path(__file__).resolve().parent
FOCUS = HERE.parent/'deadline_focus_20261006'
INVENTORY = HERE.parent/'current_inventory_20261006'
assert not (HERE/'registry.json').exists(), 'Already applied; inspect the existing registry.'
OUT = Path('/usr/xtmp/lz280/markstar_queue_alignment_20261007')/os.environ['SLURM_JOB_ID']
OUT.mkdir(parents=True, exist_ok=False)
targets = {'2rl0_flex_p4':'12827033_3','3cal_flex_p5':'12827034_4','4u3s_flex_p1':'12826673_0'}
child, old_parent, new_parent = '12827035_1','12826673_0','12812439_1'
commands = []

def run(*args):
    result = subprocess.run(args, text=True, capture_output=True, timeout=40)
    commands.append(dict(arguments=list(args), returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
    (OUT/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    result.check_returncode()
    return result.stdout

def detail(job):
    raw = run('scontrol','show','job',job,'-o')
    return dict(re.findall(r'(?:^|\s)([A-Za-z][A-Za-z0-9_]*)=(\S+)',raw))

def pending(job):
    fields = detail(job)
    assert fields['JobState']=='PENDING' and fields['RunTime']=='00:00:00' and fields.get('Restarts')=='0', fields
    assert fields['Account']=='grisman' and fields['ArrayTaskId'].isdigit(), fields
    return fields

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows):
    with path.open('w', newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')

prefs = json.loads((INVENTORY/'reporting_preferences.json').read_text())
assert set(targets)==set(prefs['omit_from_plan'])
before_raw=run('squeue','-h','-r','-u','lz280','-o','%i|%j|%T|%M|%R|%P|%E|%C')
(OUT/'queue_before.txt').write_text(before_raw)
before_queue = {r[0]:r[1:] for line in before_raw.splitlines() if (r:=line.split('|'))}
protected = {job:detail(job) for job,r in before_queue.items() if r[1]=='RUNNING' and r[0].startswith('markstar')}
for job,fields in before_queue.items():
    dependency=fields[5]
    for target in targets.values():
        if re.search(r'(?<!\d)'+re.escape(target)+r'(?!\d)',dependency):
            assert job==child and target==old_parent, (job,dependency)
before = {job:pending(job) for job in targets.values()}
before_child=pending(child)
assert 'afterany:'+old_parent in before_child['Dependency'] and before_child['Reason']!='JobHeldUser'
assert 'afterany:'+new_parent in before[old_parent]['Dependency']
assert detail(new_parent)['JobState']=='RUNNING'
save(OUT/'before.json',dict(targets=before,child=before_child,protected=protected))
paths=[FOCUS/p for p in ('design_actions.tsv','local_lane_dependencies.tsv','cpu_continuation_submissions.tsv','policy.json','registry.json')]
paths += [INVENTORY/'excluded_designs.json',INVENTORY/'reporting_preferences.json']
for path in paths:
    dest=OUT/'before'/path.parent.name/path.name
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(path,dest)

# Hold the successor while bypassing p1, before cancelling its old predecessor.
run('scontrol','hold',child)
try:
    pending(child)
    pending(old_parent)
    run('scontrol','update','JobId='+child,'Dependency=afterany:'+new_parent)
    assert 'afterany:'+new_parent in detail(child)['Dependency']
except BaseException:
    if detail(child)['JobState']=='PENDING':
        run('scontrol','update','JobId='+child,'Dependency=afterany:'+old_parent)
        run('scontrol','release',child)
    raise
actions=[]
for design,job in targets.items():
    pending(job)
    run('scancel','--state=PENDING',job)
    actions.append(dict(design=design,job=job,action='CANCELLED_NEVER_STARTED'))
    save(OUT/'actions.json',actions)
run('scontrol','release',child)
after_raw=run('squeue','-h','-r','-u','lz280','-o','%i|%j|%T|%M|%R|%P|%E|%C')
(OUT/'queue_after.txt').write_text(after_raw)
after_queue={r[0]:r[1:] for line in after_raw.splitlines() if (r:=line.split('|'))}
assert not set(targets.values()) & set(after_queue)
assert child in after_queue and 'afterany:'+new_parent in after_queue[child][5]
assert after_queue[child][3]!='(JobHeldUser)'
for job,old in protected.items():
    new=detail(job)
    assert all(new[k]==old[k] for k in ('StartTime','TimeLimit','NumCPUs','NodeList','Command')),job

definitions=read(FOCUS/'design_actions.tsv')
for row in definitions:
    if row['design'] in targets:
        row.update(action='DEFER_REMOVED_FROM_PLAN',current_status='LOCAL_PENDING_COPY_CANCELLED',current_local_job='')
write(FOCUS/'design_actions.tsv',definitions)
edges=read(FOCUS/'local_lane_dependencies.tsv')
for row in edges:
    if row['job']==child:
        row.update(parent=new_parent,action='BYPASS_OMITTED_P1_20261007')
    elif row['job'] in targets.values():
        row['action']='CANCELLED_REMOVED_FROM_PLAN'
write(FOCUS/'local_lane_dependencies.tsv',edges)
submissions=read(FOCUS/'cpu_continuation_submissions.tsv')
for row in submissions:
    if row['job']==child:
        row['after_job']=new_parent
write(FOCUS/'cpu_continuation_submissions.tsv',submissions)
exclusions=json.loads((INVENTORY/'excluded_designs.json').read_text())
exclusions['designs']=sorted(set(exclusions['designs'])|set(targets))
exclusions['local_cancelled_jobs'].update(targets)
exclusions['queue_alignment_audit']=str(OUT)
save(INVENTORY/'excluded_designs.json',exclusions)
policy=json.loads((FOCUS/'policy.json').read_text())
for field in ('primary_targets','coverage_pending','coverage_new','conditional_reserve_existing'):
    for name in targets:
        policy.get(field,{}).pop(name,None)
policy['queue_alignment_registry']=str(HERE/'registry.json')
save(FOCUS/'policy.json',policy)
prefs['scope']='Approved unified plan; omitted and conditional tiers are absent from the local queue. Unsubmitted additions remain explicitly labelled.'
save(INVENTORY/'reporting_preferences.json',prefs)

# Audit the actual local main-campaign queue against the agreed unified table.
plan=read(INVENTORY/'latest_report/unified_plan.tsv')
job_to_design={r['job']:r['design'] for r in plan if r.get('job')}
unknown=[]
active=[]
for job,r in after_queue.items():
    if not r[0].startswith('markstar') or r[6]!='64':
        continue
    if job not in job_to_design:
        unknown.append(job)
    else:
        active.append(dict(job=job,design=job_to_design[job],state=r[1],dependency=r[5]))
assert not unknown, ('Local queued tasks not present in unified plan',unknown)
record=dict(applied=datetime.now(ZoneInfo('America/New_York')).isoformat(),management_job=os.environ['SLURM_JOB_ID'],
            audit=str(OUT),cancelled=actions,dependency_change=dict(job=child,old_parent=old_parent,new_parent=new_parent),
            protected_running=len(protected),active_local_main=active,local_queue_outside_plan=[],
            not_submitted=[r['design'] for r in plan if r['state']=='待提交'],
            dcc_status='HANDOFF_ONLY_NO_REMOTE_ACCESS_OR_CURRENT_JOB_MAPPING',
            scope='All actual local main-campaign jobs are retained plan entries; future unsubmitted plan entries remain labelled.')
save(HERE/'registry.json',record)
old_registry=json.loads((FOCUS/'registry.json').read_text())
old_registry['queue_alignment_registry']=str(HERE/'registry.json')
save(FOCUS/'registry.json',old_registry)
print(json.dumps(record,ensure_ascii=False,indent=2),flush=True)
