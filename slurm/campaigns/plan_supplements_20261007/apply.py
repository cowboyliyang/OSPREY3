"""Apply authorized plan additions, pending-only deferral, and one GPU submission."""
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime
from zoneinfo import ZoneInfo

assert os.environ.get('SLURM_JOB_ID')
HERE = Path(__file__).resolve().parent
FOCUS = HERE.parent/'deadline_focus_20261006'
INVENTORY = HERE.parent/'current_inventory_20261006'
OUT = Path('/usr/xtmp/lz280/markstar_plan_supplements_20261007')/('apply_'+os.environ['SLURM_JOB_ID'])
OUT.mkdir(parents=True, exist_ok=False)
registry = json.loads((HERE/'registry.json').read_text())
prep = Path(registry['prep_root'])
assert (prep/'READY').exists()
assert not registry.get('pack_job'), 'Do not submit a duplicate GPU job.'
targets = {'3eb6_flex_p7':'12827030_0', '4wem_flex_p12':'12827029_0', '5dc4_flex_p6':'12826673_3'}

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')

def queue():
    raw = subprocess.check_output(['squeue','-h','-r','-u','lz280','-o','%i|%T|%E'], text=True, timeout=40)
    return raw, {r[0]:r[1:] for line in raw.splitlines() if len(r:=line.split('|'))==3}

changed = [FOCUS/p for p in ('dcc_new_designs.tsv','dcc_new_cases.txt','dcc_launch_priority.tsv','dcc_actions.tsv','design_actions.tsv','policy.json','registry.json')]
changed += [INVENTORY/'reporting_preferences.json']
for path in changed:
    backup = OUT/'before'/path.parent.name/path.name
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup)
raw, before = queue()
(OUT/'queue_before.txt').write_text(raw)
for job, (_, dependency) in before.items():
    assert not any(re.search(r'(?<!\d)'+re.escape(target)+r'(?!\d)', dependency) for target in targets.values()), (job,dependency)
actions = []
for design, job in targets.items():
    query = subprocess.run(['scontrol','show','job',job,'-o'], text=True, capture_output=True, check=True)
    (OUT/(job+'_before.txt')).write_text(query.stdout)
    info = dict(re.findall(r'(\w+)=([^\s]+)', query.stdout))
    action = 'PRESERVED_STARTED_OR_FINISHED'
    if info.get('JobState')=='PENDING' and info.get('RunTime')=='00:00:00' and info.get('Restarts')=='0':
        assert info.get('ArrayTaskId','').isdigit()
        subprocess.run(['scancel','--state=PENDING',job], check=True)
        action = 'CANCEL_REQUESTED_IF_STILL_PENDING'
    actions.append(dict(design=design, job=job, state_before=info.get('JobState'), action=action,
                        prior_dependency=info.get('Dependency','')))
    save(OUT/'conditional_actions.json', actions)
raw, after = queue()
(OUT/'queue_after.txt').write_text(raw)
deferred = {}
for row in actions:
    live = after.get(row['job'])
    if live is None and row['action']=='CANCEL_REQUESTED_IF_STILL_PENDING':
        row['action']='CANCELLED_PENDING'
        deferred[row['design']] = dict(previous_job=row['job'], state='条件保留，暂不排队', audit=str(OUT))
    elif live and live[0]=='RUNNING':
        row['action']='PRESERVED_STARTED_DURING_UPDATE'
    else:
        raise RuntimeError(('Unexpected conditional job state',row,live))
save(OUT/'conditional_actions.json', actions)
save(INVENTORY/'conditional_deferred.json', dict(updated=datetime.now(ZoneInfo('America/New_York')).isoformat(),
     designs=deferred, automatic_resubmission=False, preserve_started=True))

new = next(r for r in read(HERE/'designs.tsv') if r['design_id']=='5dc0_flex_p9')
dcc = read(FOCUS/'dcc_new_designs.tsv')
assert not any(r['design_id']==new['design_id'] for r in dcc)
anchor = next(r for r in dcc if r['design_id']=='5dc0_flex_p10')
handoff = dict(anchor)
handoff.update(new)
handoff.update(group_index='0', source_package=str(prep/'package'),
               handoff_index=str(max(int(r['handoff_index']) for r in dcc)+1))
dcc.append(handoff)
write(FOCUS/'dcc_new_designs.tsv', dcc)
(FOCUS/'dcc_new_cases.txt').write_text(''.join(r['design_id']+'\n' for r in dcc))
priority = read(FOCUS/'dcc_launch_priority.tsv')
priority.append(dict(priority=str(max(int(r['priority']) for r in priority)+1), design_id=new['design_id'],
                     role='USER_ADDITION', prior_submission='NEW_DEFINITION_NOT_SUBMITTED_REMOTELY', gpu_dependency='NONE',
                     runtime_evidence='User selected p9; both input preflights passed; MARK and PACK times unmeasured'))
write(FOCUS/'dcc_launch_priority.tsv', priority)
dcc_actions = read(FOCUS/'dcc_actions.tsv')
dcc_actions.append(dict(design_id=new['design_id'], action='KEEP_USER_ADDITION', reported_status='NOT_SUBMITTED_REMOTELY',
                        requested_cpus='64', limit_hours='336', combined_concurrency='6', gate='NONE',
                        cpu_rule='Match 5dc0_p4: Intel Xeon Platinum 8462Y+; verify actual model before launch'))
write(FOCUS/'dcc_actions.tsv', dcc_actions)
design_actions = read(FOCUS/'design_actions.tsv')
for row in design_actions:
    if row['design'] in deferred:
        row.update(action='CONDITIONAL_PLAN_ONLY', current_status='NOT_QUEUED', current_local_job='')
design_actions.append(dict(design=new['design_id'], system='5dc0', action='KEEP_USER_ADDITION', current_assignment='dcc',
                           proposed_destination='dcc', current_local_job='', current_status='NOT_SUBMITTED_REMOTELY',
                           old_mark_hours='', old_pack_minutes='', historical_ratio='NA'))
for definition in read(HERE/'designs.tsv'):
    if definition['design_id']==new['design_id']:
        continue
    assert not any(r['design']==definition['design_id'] for r in design_actions)
    design_actions.append(dict(design=definition['design_id'], system=definition['system'], action='PLAN_ADDITION_READY',
                               current_assignment='local', proposed_destination='compsci', current_local_job='',
                               current_status='NOT_SUBMITTED', old_mark_hours='', old_pack_minutes='', historical_ratio='NA'))
write(FOCUS/'design_actions.tsv', design_actions)
policy = json.loads((FOCUS/'policy.json').read_text())
for name in deferred:
    policy['primary_targets'].pop(name, None)
policy['conditional_plan_only'] = deferred
policy.setdefault('supplemental_targets',{})[new['design_id']] = 'dcc'
save(FOCUS/'policy.json', policy)
prefs = json.loads((INVENTORY/'reporting_preferences.json').read_text())
for addition in prefs['plan_additions']:
    if addition['design'] in ('3k3q_flex_m1','5d68_flex_p4'):
        addition['definition_exists']=True
        addition['source'] += '; MARK and PACK input preflights passed in job '+registry['prep_job']
prefs['plan_additions'].append(dict(design=new['design_id'], definition_exists=True, destination='dcc', source='User explicitly requested DCC p9 on October 7'))
prefs['focus_additions']='5dc0 p9 added to DCC handoff; retain p4 and p10; no remote submission receipt yet.'
prefs['conditional_scheduling']='Table annotation only; never-started submissions cancelled, no automatic resubmission.'
save(INVENTORY/'reporting_preferences.json', prefs)
receipt = subprocess.check_output(['sbatch','--parsable','--account=grisman',
                '--export=ALL,EXTENSION_PREP='+str(prep), str(prep/'run_group.slurm')], text=True).strip()
job = receipt.split(';')[0]
assert job.isdigit(), receipt
save(OUT/'pack_submission.json', dict(job=job, receipt=receipt, script=str(prep/'run_group.slurm')))
registry.update(pack_job=job, packstar_status='SUBMITTED', pack_output='/usr/xtmp/lz280/markstar_plan_supplements_20261007/pack_'+job,
                apply_job=os.environ['SLURM_JOB_ID'], apply_audit=str(OUT), conditional_actions=actions,
                dcc_addition=dict(design='5dc0_flex_p9', status='HANDOFF_READY_NO_REMOTE_RECEIPT', requested_cpus=64,
                                  limit_hours=336, combined_concurrency=6, cpu_model='Intel Xeon Platinum 8462Y+'))
save(HERE/'registry.json', registry)
focus_registry = json.loads((FOCUS/'registry.json').read_text())
focus_registry['october7_supplements']=dict(registry=str(HERE/'registry.json'), dcc_addition=registry['dcc_addition'], conditional_actions=actions)
save(FOCUS/'registry.json', focus_registry)
print(json.dumps(registry, ensure_ascii=False, indent=2), flush=True)
