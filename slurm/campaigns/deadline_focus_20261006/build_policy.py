"""Write the active focus registry and a pending-only DCC handoff, through Slurm."""
import collections
import csv
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here=Path(__file__).parent
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rs):
    with p.open('w') as f:
        fields=list(rs[0])
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rs:
            row=dict(r)
            if row[fields[-1]] in ('',None):row[fields[-1]]='NA'
            w.writerow(row)
rows=read(Path('/usr/xtmp/lz280/markstar_deadline_proposal_12826811/design_actions.tsv'))
registry=json.loads((here/'registry.json').read_text())
registry.update(packstar_status='SUBMITTED',markstar_status='TWO_NEW_LOCAL_JOBS_SUBMITTED_WITH_PACKSTAR_GATE',
                local_lane_limits=dict(compsci=8,fennario=4),dcc_concurrency=6)
(here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
local={r['design']:r for r in read(here/'local_actions.tsv')}
newlocal={r['design']:r for r in read(here/'new_local_submissions.tsv')}
remote=[];allrows=[]
for r in rows:
    r=dict(r)
    if r['design']=='2xxm_flex_p3':r['action']='HOLD_COVERAGE_FALLBACK'
    if r['design'] in newlocal:
        r['current_local_job']=newlocal[r['design']]['job'];r['current_status']='SUBMITTED_WITH_PACKSTAR_GATE'
    elif r['design'] in local:
        actual=local[r['design']]
        r['current_status']=actual['after']
        if actual['action']=='PROTECT_STARTED':r['action']='KEEP_RUNNING'
        if actual['action'].startswith('RETIRE'):r['current_status']='LOCAL_PENDING_COPY_CANCELLED'
    allrows.append(r)
    if r['current_assignment']=='dcc' or r['proposed_destination']=='dcc':
        d=dict(design_id=r['design'],action=r['action'],reported_status=r['current_status'],
               requested_cpus=64,limit_hours=336,combined_concurrency=6,
               gate='PACKSTAR_NORMAL_COMPLETION_LE_1800_SECONDS' if r['action']=='PROPOSE_NEW_PRIMARY' else '',
               cpu_rule='Use same actual CPU model as existing current-round same-system anchor; verify model before launch')
        remote.append(d)
write(here/'design_actions.tsv',allrows)
write(here/'dcc_actions.tsv',remote)
for filename,predicate in [
 ('dcc_cancel_if_pending_cases.txt',lambda r:r['action'].startswith('DEFER')),
 ('dcc_keep_existing_pending_cases.txt',lambda r:r['action'] in ('KEEP_COVERAGE','KEEP_PRIMARY')),
 ('dcc_protected_cases.txt',lambda r:r['action'] in ('KEEP_COMPLETED','KEEP_REPORTED_RUNNING','KEEP_PREVIOUSLY_STARTED_STATUS_UNCONFIRMED'))]:
    (here/filename).write_text(''.join(r['design_id']+'\n' for r in remote if predicate(r)))
inv=json.loads(Path('/usr/xtmp/lz280/markstar_deadline_coverage_12826781/inventory.json').read_text())
coverage=read(Path('/usr/xtmp/lz280/markstar_deadline_proposal_12826811/system_coverage.tsv'))
write(here/'system_coverage.tsv',coverage)
(here/'external_status_snapshot.json').write_text(json.dumps(inv['external_status'],indent=2)+'\n')
cfg=json.loads(Path('/home/users/lz280/notes/markstar_deadline_proposal_20261006.json').read_text())
cfg.update(status='LOCAL_APPLIED_DCC_HANDOFF_NOT_YET_APPLIED',user_scope='All unstarted work can be reprioritized; all started work and completed results remain protected.',
           snapshot_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
           dcc_counts=dict(collections.Counter(r['action'] for r in remote)),
           local_adjustment=json.loads((here/'local_adjustment.json').read_text()),
           local_submission_registry=str(here/'registry.json'))
(here/'policy.json').write_text(json.dumps(cfg,indent=2)+'\n')
old_registry_path=here.parent/'frontier_extension7_20261006/registry.json'
old_registry=json.loads(old_registry_path.read_text())
old_registry.update(superseded_by='../deadline_focus_20261006/README.md',
    current_local_actions='../deadline_focus_20261006/local_actions.tsv',
    current_dcc_actions='../deadline_focus_20261006/dcc_actions.tsv',
    packstar_status='PENDING_JOB_REPLACED',packstar_replacement_job=registry['packstar_job'],
    markstar_status='PENDING_QUEUE_PRUNED_SEE_DEADLINE_FOCUS')
old_registry_path.write_text(json.dumps(old_registry,indent=2)+'\n')
old_handoff_path=here.parent/'frontier_extension7_20261006/dcc_handoff.json'
old_handoff=json.loads(old_handoff_path.read_text())
old_handoff.update(status='SUBMISSION_CONFIRMED_BY_USER_SUBSEQUENT_PENDING_REPLAN',
    superseded_by='../deadline_focus_20261006/README.md',
    held_local_copy=dict(design='4wwi_flex_p1',job='12812439_24',status='CANCELLED_PENDING_DUPLICATE_AFTER_ACCEPTANCE'))
old_handoff_path.write_text(json.dumps(old_handoff,indent=2)+'\n')
assert len(remote)==104 and len(allrows)==169
assert len(coverage)==38
assert all(s['completed'] or s['running'] or s['pending_or_new'] for s in coverage)
for path in here.glob('*.py'):compile(path.read_text(),str(path),'exec')
for path in here.glob('*.slurm'):subprocess.run(['bash','-n',str(path)],check=True)
subprocess.run(['git','diff','--check'],cwd=here.parents[2],check=True)
print('DCC_ACTION_COUNTS',dict(collections.Counter(r['action'] for r in remote)))
print('LOCAL_RUNNING_PROTECTED',12,'EXISTING_PENDING_KEEP',11,'NEW_GATED_LOCAL',2,'HELD_FALLBACK',1)
print('ALL_38_SYSTEMS_HAVE_COVERAGE_OR_RETAINED_CANDIDATE')
