"""Remove the GPU gate and queue useful preflighted continuation designs.

Only pending jobs are replaced. Existing frozen packages and running jobs
remain untouched. Slurm may still schedule other users between allocations.
"""
import collections
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here=Path(__file__).parent
registry=json.loads((here/'registry.json').read_text())
assert 'cpu_continuation' not in registry, 'Already applied; inspect registry instead of resubmitting'
out=Path('/usr/xtmp/lz280/markstar_deadline_focus_20261006')/('cpu_continuation_'+os.environ['SLURM_JOB_ID'])
package=out/'package';package.mkdir(parents=True,exist_ok=False)
sources=[Path(registry['markstar_prep']),Path('/usr/xtmp/lz280/frontier_extension7_20261006/markstar_prep_12826665')]
available={}
for source in sources:
    assert (source/'READY').exists(),source
    available.update(json.loads((source/'package/run_configs.json').read_text()))
requests=[
 ('4wem_flex_p12','fennario','12812439_2','12826957_0','PRIMARY_NO_GPU_GATE'),
 ('3eb6_flex_p7','compsci','12812439_26','12826958_0','PRIMARY_NO_GPU_GATE'),
 ('1a0r_flex_p0','compsci','12812439_52','12826674_6','CONTINUATION_RESERVE'),
 ('4wyq_flex_p3','compsci','12812439_41','12826674_4','CONTINUATION_RESERVE'),
 ('2rl0_flex_p4','compsci','12826674_1','12826674_2','CONTINUATION_RESERVE'),
 ('3cal_flex_p5','compsci','12826674_7','12826674_8','CONTINUATION_RESERVE'),
 ('4u3s_flex_p2','fennario','12826673_0','12826673_2','CONTINUATION_RESERVE')]

def query():
    raw=subprocess.check_output(['squeue','-u','lz280','-r','-h','-o','%i|%T|%E|%R'],text=True)
    return raw,{p[0]:p[1:] for line in raw.splitlines() if (p:=line.split('|'))}
def write(path,rows):
    keys=list(rows[0])
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:
            row=dict(row)
            if row[keys[-1]] in ('',None):row[keys[-1]]='NA'
            w.writerow(row)
def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))

raw,before=query();(out/'queue_before.txt').write_text(raw)
protected={j for j,v in before.items() if v[0]=='RUNNING'}
assert all(old in before and before[old][0]=='PENDING' for _,_,_,old,_ in requests[:2])
assert all(old not in before for _,_,_,old,_ in requests[2:]),'A reserve already has an active copy'
configs={};rows=[];indices={'fennario':0,'compsci':0}
for name,group,parent,old,role in requests:
    config=json.loads(json.dumps(available[name]));row=config['row']
    assert row['destination']==group
    row['group_index']=indices[group];indices[group]+=1
    assert int(row['requested_cpus'])==64 and int(row['java_heap_gib'])<192
    assert hashlib.sha256(Path(config['input_pdb']).read_bytes()).hexdigest()==row['pdb_sha256']
    assert '-Dosprey.bench.numCPUs=64' in config['command_template']
    assert '-XX:ActiveProcessorCount=64' in config['command_template']
    assert all(Path(p).exists() for p in config['command_template'][config['command_template'].index('-cp')+1].split(':'))
    configs[name]=config
    rows.append({k:row[k] for k in ['design_id','kind','destination','group_index','expected_sequences','pdb_sha256','java_heap_gib','reference_design','reference_mark_path']})
write(package/'designs.tsv',rows)
(package/'run_configs.json').write_text(json.dumps(configs,indent=2)+'\n')
for filename in ['markstar_runner.py','markstar_array.slurm']:
    source=here/filename
    if filename.endswith('.py'):
        text=source.read_text();compile(text,str(source),'exec')
        assert 'FOCUS_PACK_ROOT' not in text and 'screening.json' not in text
    else:subprocess.run(['bash','-n',str(source)],check=True)
    shutil.copy2(source,out/filename)
(out/'READY').write_text('Seven preserved preflighted commands; CPU execution independent of GPU results\n')
submitted=[];held_old=[];released=False
try:
    for _,_,_,old,_ in requests[:2]:
        _,live=query();assert live[old][0]=='PENDING',old
        subprocess.run(['scontrol','hold',old],check=True);held_old.append(old)
        detail=subprocess.check_output(['scontrol','show','job',old,'-o'],text=True)
        assert 'JobState=PENDING ' in detail and 'RunTime=00:00:00 ' in detail and 'Restarts=0 ' in detail
        (out/(old+'_before.txt')).write_text(detail)
    for name,group,parent,old,role in requests:
        index=configs[name]['row']['group_index']
        args=['sbatch','--parsable','--hold','--account=grisman','--job-name=markstar-cont-'+group,
              '--partition='+('grisman' if group=='fennario' else 'compsci'),'--array='+str(index)+'-'+str(index)+'%1',
              '--dependency=afterany:'+parent,
              '--export=ALL,MARKSTAR_EXTENSION_PREP='+str(out)+',MARKSTAR_EXTENSION_GROUP='+group]
        if group=='fennario':args+=['--constraint=a5000','--exclude=grisman-37,grisman-40,jerry[1-7]']
        else:args+=['--exclude=linux[31-40]']
        args.append(str(out/'markstar_array.slurm'))
        jid=subprocess.check_output(args,text=True).strip().split(';')[0];assert jid.isdigit()
        record=dict(design=name,job=jid+'_'+str(index),group=group,after_job=parent,
                    replaced_or_previous_job=old,role=role,gpu_dependency='NONE',prep=str(out))
        submitted.append(record)
        (out/'submitted.json').write_text(json.dumps(submitted,indent=2)+'\n')
        detail=subprocess.check_output(['scontrol','show','job',record['job'],'-o'],text=True)
        assert 'Account=grisman ' in detail and 'CPUs/Task=64 ' in detail and 'TimeLimit=14-00:00:00 ' in detail
        assert '12826950' not in detail and 'JobHeldUser' in detail
    for old in held_old:
        _,live=query();assert live[old][0]=='PENDING'
    for old in held_old:subprocess.run(['scancel','--state=PENDING',old],check=True)
    for record in submitted:subprocess.run(['scontrol','release',record['job']],check=True)
    released=True
except BaseException:
    if not released:
        for record in submitted:subprocess.run(['scancel','--state=PENDING',record['job']],check=False)
        for old in held_old:subprocess.run(['scontrol','release',old],check=False)
    raise
raw,after=query();(out/'queue_after.txt').write_text(raw)
for r in submitted:
    assert r['job'] in after and '12826950' not in after[r['job']][1]
assert all(old not in after for old in held_old)
write(here/'cpu_continuation_submissions.tsv',submitted)
write(here/'new_local_submissions.tsv',submitted[:2])
local_rows=read(here/'local_lane_dependencies.tsv')
for r in submitted:local_rows.append(dict(job=r['job'],parent=r['after_job'],action='CPU_ONLY_CONTINUATION'))
write(here/'local_lane_dependencies.tsv',local_rows)
remote_reserves={'4wyu_flex_p1','3bua_flex_m2','3ma2_flex_p1','4wwi_flex_p2'}
remote=read(here/'dcc_actions.tsv')
for r in remote:
    r['gate']='NONE'
    if r['design_id'] in remote_reserves:r['action']='KEEP_CONTINUATION_RESERVE'
write(here/'dcc_actions.tsv',remote)
(here/'dcc_cancel_if_pending_cases.txt').write_text(''.join(r['design_id']+'\n' for r in remote if r['action'].startswith('DEFER')))
(here/'dcc_keep_existing_pending_cases.txt').write_text(''.join(r['design_id']+'\n' for r in remote if r['action'] in ('KEEP_COVERAGE','KEEP_PRIMARY','KEEP_CONTINUATION_RESERVE')))
actions=read(here/'design_actions.tsv');bydesign={r['design']:r for r in submitted}
for r in actions:
    if r['design'] in bydesign:
        job=bydesign[r['design']]
        r.update(current_local_job=job['job'],current_status='SUBMITTED_CPU_ONLY',
                 action='KEEP_CONTINUATION_RESERVE' if job['role']=='CONTINUATION_RESERVE' else 'KEEP_PRIMARY')
    elif r['design'] in remote_reserves:r['action']='KEEP_CONTINUATION_RESERVE'
write(here/'design_actions.tsv',actions)
record=dict(applied_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audit=str(out),
            gpu_dependency=None,local_submitted=submitted,retired_gpu_gated_jobs=held_old,
            added_local_reserves=5,added_dcc_reserves=sorted(remote_reserves),dcc_retained_existing_pending=10,
            dcc_cancel_if_pending=57,dcc_new=4,protected_running_snapshot=sorted(protected),
            scheduler_limitation='Dependencies prepare successors; Slurm cannot guarantee uninterrupted ownership across separate allocations.')
registry.update(markstar_status='CPU_ONLY_NO_GPU_GATE',markstar_prep=str(out),new_local_jobs=submitted[:2],
                cpu_continuation=record,dcc_status='UPDATED_HANDOFF_NOT_APPLIED_REMOTELY')
(here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
policy=json.loads((here/'policy.json').read_text())
policy.update(status='CPU_CONTINUATION_LOCAL_APPLIED_DCC_HANDOFF_PENDING',gpu_precondition=False,
              cpu_continuation=record,conditional_reserve_existing={},
              dcc_counts=dict(collections.Counter(r['action'] for r in remote)),
              rules=[r for r in policy['rules'] if not any(word in r for word in ['Measure new PACK','Prefer PACK','One primary','After October 16'])]+[
                  'User prioritizes continuous useful CPU work; launch preflighted MARK* independently of GPU timings.',
                  'PACK* timings inform interpretation and later tier choice, never block these queued CPU tasks.',
                  'Retain compsci/fennario/DCC total limits of 8/4/6; preserve all started work.',
                  'Keep high-potential continuation tiers queued after CPU predecessors; do not restore the low-value intermediate backlog.',
                  'October 16 remains a planning cutoff for a full 14-day observation by October 30, not a GPU-dependent runtime assertion.'])
(here/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
assert len((here/'dcc_cancel_if_pending_cases.txt').read_text().splitlines())==57
assert len((here/'dcc_keep_existing_pending_cases.txt').read_text().splitlines())==10
subprocess.run(['git','diff','--check'],cwd=here.parents[2],check=True)
(out/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2),flush=True)
