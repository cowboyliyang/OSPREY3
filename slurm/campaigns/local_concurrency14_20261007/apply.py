"""Raise local MARK* concurrency to 9+5 by splitting two serial task chains."""
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

assert os.environ.get('SLURM_JOB_ID'), 'Submit with account grisman.'
HERE=Path(__file__).resolve().parent
FOCUS=HERE.parent/'deadline_focus_20261006'
INVENTORY=HERE.parent/'current_inventory_20261006'
PLAN=json.loads((HERE/'plan.json').read_text())
assert not (HERE/'registry.json').exists(), 'Inspect the existing application instead of repeating it.'
OUT=Path('/usr/xtmp/lz280/markstar_local_concurrency14_20261007')/os.environ['SLURM_JOB_ID']
OUT.mkdir(parents=True,exist_ok=False)
commands=[]

def save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def run(*args):
    result=subprocess.run(args,capture_output=True,text=True,timeout=40)
    commands.append(dict(arguments=list(args),returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    save(OUT/'commands.json',commands)
    result.check_returncode()
    return result.stdout

def state(job):
    return dict(re.findall(r'(?:^|\s)([A-Za-z][A-Za-z0-9_]*)=(\S+)',run('scontrol','show','job',job,'-o')))

def pending(change):
    fields=state(change['job'])
    assert fields['JobState']=='PENDING' and fields['RunTime']=='00:00:00' and fields.get('Restarts')=='0',fields
    assert fields['Account']=='grisman' and fields['NumCPUs']=='64' and fields['Partition']==change['partition'],fields
    assert fields.get('ArrayTaskId','').isdigit(),fields
    return fields

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))

def write(path,rows):
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def queue():
    raw=run('squeue','-h','-r','-u','lz280','-o','%i|%j|%T|%P|%E|%C|%a')
    rows={}
    for line in raw.splitlines():
        p=line.split('|')
        if len(p)==7 and p[1].startswith('markstar') and p[5]=='64':
            rows[p[0]]=dict(zip(('name','state','partition','dependency','cpus','account'),p[1:]))
    return raw,rows

def parents(rows):
    links={}
    for job,row in rows.items():
        if row['dependency'] in ('(null)','', '0'):
            links[job]=None
            continue
        match=re.fullmatch(r'afterany:(\d+(?:_\d+)?)(?:\(unfulfilled\))?',row['dependency'])
        assert match,(job,row)
        links[job]=match[1] if match[1] in rows else None
    return links

def validate(rows,links,limits):
    children=Counter(p for p in links.values() if p is not None)
    assert all(n<=1 for n in children.values()), 'A serial lane must not branch.'
    roots=Counter(rows[j]['partition'] for j,p in links.items() if p is None)
    assert roots['compsci']<=limits['compsci'] and roots['grisman']<=limits['fennario'],roots
    for job,parent in links.items():
        if parent:
            assert rows[job]['partition']==rows[parent]['partition']
        seen=set()
        cursor=job
        while cursor is not None:
            assert cursor not in seen,'Dependency cycle.'
            seen.add(cursor)
            cursor=links[cursor]
    return dict(roots)

raw,before_queue=queue()
(OUT/'queue_before.txt').write_text(raw)
before_links=parents(before_queue)
before_roots=validate(before_queue,before_links,PLAN['old_limits'])
table={r['design']:r for r in read(INVENTORY/'current_plan.tsv')}
prefs=json.loads((INVENTORY/'reporting_preferences.json').read_text())
excluded=set(json.loads((INVENTORY/'excluded_designs.json').read_text())['designs'])
before={}
for change in PLAN['changes']:
    assert table[change['design']]['job']==change['job'] and table[change['design']]['recommendation']=='保留'
    assert change['design'] not in set(prefs['conditional_designs'])|excluded
    before[change['job']]=pending(change)
    assert before[change['job']]['Reason']!='JobHeldUser'
    assert before_links[change['job']]==change['old_parent']
planned_links=dict(before_links)
for change in PLAN['changes']:
    planned_links[change['job']]=None
planned_roots=validate(before_queue,planned_links,PLAN['new_limits'])
assert planned_roots=={'compsci':9,'grisman':5},planned_roots
protected={job:state(job) for job,row in before_queue.items() if row['state']=='RUNNING'}
save(OUT/'before.json',dict(selected=before,running=protected,lane_roots=before_roots))
for path in (FOCUS/'policy.json',FOCUS/'registry.json',FOCUS/'local_lane_dependencies.tsv',
             FOCUS/'cpu_continuation_submissions.tsv',FOCUS/'design_actions.tsv',INVENTORY/'reporting_preferences.json'):
    backup=OUT/'before'/path.parent.name/path.name
    backup.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(path,backup)

held=[]
updated=[]
try:
    for change in PLAN['changes']:
        pending(change)
        run('scontrol','hold',change['job'])
        held.append(change)
        pending(change)
    for change in PLAN['changes']:
        pending(change)
        run('scontrol','update','JobId='+change['job'],'Dependency=')
        updated.append(change)
        assert state(change['job'])['Dependency'] in ('(null)','', '0')
except BaseException:
    rollback=True
    for change in reversed(updated):
        try:
            pending(change)
            run('scontrol','update','JobId='+change['job'],'Dependency=afterany:'+change['old_parent'])
        except BaseException:
            rollback=False
    if rollback:
        for change in held:
            run('scontrol','release',change['job'])
    save(OUT/'failure.json',dict(rollback_succeeded=rollback,held=held))
    raise
for change in held:
    run('scontrol','release',change['job'])
raw,after_queue=queue()
(OUT/'queue_after.txt').write_text(raw)
after_roots=validate(after_queue,parents(after_queue),PLAN['new_limits'])
for change in PLAN['changes']:
    info=state(change['job'])
    assert info['Dependency'] in ('(null)','', '0') and info['Reason']!='JobHeldUser',info
    assert all(info[k]==before[change['job']][k] for k in ('Command','NumCPUs','MinMemoryNode','TimeLimit','Partition','Account')),info
for job,old in protected.items():
    current=state(job)
    assert all(current[k]==old[k] for k in ('StartTime','Command','NumCPUs','NodeList','TimeLimit')),job

changed={r['job']:r for r in PLAN['changes']}
edges=read(FOCUS/'local_lane_dependencies.tsv')
for row in edges:
    if row['job'] in changed:
        row.update(parent='NONE',action='INDEPENDENT_LANE_LOCAL_CONCURRENCY14')
write(FOCUS/'local_lane_dependencies.tsv',edges)
submissions=read(FOCUS/'cpu_continuation_submissions.tsv')
for row in submissions:
    if row['job'] in changed:
        row['after_job']='NONE'
write(FOCUS/'cpu_continuation_submissions.tsv',submissions)
actions=read(FOCUS/'design_actions.tsv')
for row in actions:
    if row['current_local_job'] in changed:
        row['current_status']='RELEASED_INDEPENDENT_LANE'
write(FOCUS/'design_actions.tsv',actions)
policy=json.loads((FOCUS/'policy.json').read_text())
policy['slots']=PLAN['new_limits']
policy['local_concurrency_registry']=str(HERE/'registry.json')
policy['rules']=[r.replace('8/4/6','9/5/6') for r in policy['rules']]
save(FOCUS/'policy.json',policy)
registry=json.loads((FOCUS/'registry.json').read_text())
registry['local_concurrency_registry']=str(HERE/'registry.json')
save(FOCUS/'registry.json',registry)
prefs['concurrency_limits']=PLAN['new_limits']
save(INVENTORY/'reporting_preferences.json',prefs)
record=dict(PLAN,status='APPLIED',management_job=os.environ['SLURM_JOB_ID'],audit=str(OUT),
            applied=datetime.now(ZoneInfo('America/New_York')).isoformat(),
            protected_running_tasks=len(protected),lane_roots_before=before_roots,lane_roots_after=after_roots,
            scheduler_states={r['job']:state(r['job'])['JobState'] for r in PLAN['changes']},
            benchmark_jobs_submitted=[],dcc_changed=False,conditional_jobs_restarted=[])
save(HERE/'registry.json',record)
print(json.dumps(record,ensure_ascii=False,indent=2),flush=True)
