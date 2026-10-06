"""Retire only live PENDING MARK* tasks outside the approved deadline focus."""
import csv
import datetime
import json
import os
from pathlib import Path
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm'
here = Path(__file__).parent
out = Path('/usr/xtmp/lz280/markstar_deadline_focus_20261006') / ('apply_' + os.environ['SLURM_JOB_ID'])
out.mkdir(parents=True, exist_ok=False)
rows = list(csv.DictReader(Path('/usr/xtmp/lz280/markstar_deadline_proposal_12826811/design_actions.tsv').open(), delimiter='\t'))
byjob = {r['current_local_job']: r for r in rows if r['current_local_job']}
arrays = '12812439,12826673,12826674'

def query():
    raw = subprocess.check_output(['squeue', '-r', '-h', '-j', arrays, '-o', '%i|%T|%R|%P'], text=True)
    return raw, {p[0]: dict(state=p[1], reason=p[2], partition=p[3]) for line in raw.splitlines() if (p := line.split('|'))}

def write_tsv(path, records):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(records)

raw, before = query()
(out/'queue_before.txt').write_text(raw)
assert set(before) <= set(byjob), set(before) - set(byjob)
protected = {j for j, v in before.items() if v['state'] != 'PENDING'}
decisions = []
for job, state in before.items():
    row = byjob[job]
    if job in protected:
        action = 'PROTECT_STARTED'
    elif row['design'] == '4wwi_flex_p1':
        action = 'RETIRE_DUPLICATE_DCC_ACCEPTED'
    elif row['design'] == '2xxm_flex_p3':
        action = 'HOLD_COVERAGE_FALLBACK'
    elif row['action'].startswith('KEEP'):
        action = 'KEEP_PENDING'
    else:
        assert row['action'].startswith('DEFER'), row
        action = 'RETIRE_PENDING'
    decisions.append(dict(job=job, design=row['design'], action=action, before=state['state'], after=''))
write_tsv(out/'planned_actions.tsv', decisions)
for r in decisions:
    if r['action'] == 'PROTECT_STARTED':
        continue
    raw, live = query()
    if r['job'] not in live or live[r['job']]['state'] != 'PENDING':
        r['action'] = 'SKIP_STATE_CHANGED'
        continue
    details = subprocess.check_output(['scontrol', 'show', 'job', r['job'], '-o'], text=True)
    (out/(r['job']+'_before.txt')).write_text(details)
    if 'Restarts=0 ' not in details or 'RunTime=00:00:00 ' not in details:
        r['action'] = 'PROTECT_PREVIOUS_EXECUTION'
        continue
    if r['action'].startswith('RETIRE'):
        # Server-side state filtering also protects jobs starting after the query.
        subprocess.run(['scancel', '--state=PENDING', r['job']], check=True)
    elif r['action'] == 'HOLD_COVERAGE_FALLBACK':
        subprocess.run(['scontrol', 'hold', r['job']], check=True)
    elif r['action'] == 'KEEP_PENDING' and r['job'].startswith('12812439_'):
        # Only the three uncovered systems remain in the original local queue.
        subprocess.run(['scontrol', 'update', 'JobId='+r['job'], 'ExcNodeList=linux[31-40]'], check=True)
    print(r['action'], r['design'], r['job'], flush=True)
raw, after = query()
(out/'queue_after.txt').write_text(raw)
for r in decisions:
    r['after'] = after.get(r['job'], {}).get('state', 'ABSENT_FROM_ACTIVE_QUEUE')
    if r['job'] in protected:
        assert r['after'] != 'PENDING'
    if r['action'].startswith('RETIRE'):
        assert r['after'] != 'PENDING', r
write_tsv(out/'local_actions.tsv', decisions)
write_tsv(here/'local_actions.tsv', decisions)
summary = dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               audit=str(out), user_scope='Only started work is protected; unstarted work on all clusters is reprioritized.',
               running_before=len(protected), pending_before=sum(v['state']=='PENDING' for v in before.values()),
               retired_pending=sum(r['action'].startswith('RETIRE') and r['after']=='ABSENT_FROM_ACTIVE_QUEUE' for r in decisions),
               kept_pending=[r for r in decisions if r['action']=='KEEP_PENDING'],
               held_fallback=[r for r in decisions if r['action']=='HOLD_COVERAGE_FALLBACK'],
               new_tiers_status='NOT_YET_PREPARED', remote_status='HANDOFF_REQUIRED_NO_REMOTE_MUTATION')
(here/'local_adjustment.json').write_text(json.dumps(summary, indent=2)+'\n')
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2), flush=True)
