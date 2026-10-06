"""Submit the approved local arrays and hold one pending DCC handoff copy."""
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import time

assert os.environ.get('SLURM_JOB_ID')
repo=Path('/home/users/lz280/OSPREY3-fresh-packstar')
here=repo/'slurm/campaigns/frontier_extension7_20261006'
campaign=repo/'slurm/campaigns/markstar_dcc_20261005'
prep=Path(os.environ['MARKSTAR_EXTENSION_PREP'])
assert (prep/'READY').is_file()
out=Path('/usr/xtmp/lz280/frontier_extension7_20261006')/('markstar_launch_'+os.environ['SLURM_JOB_ID'])
out.mkdir()
registry=json.loads((here/'registry.json').read_text())
assert registry['markstar_status']=='PLAN_ONLY','Refuse duplicate submissions'
def read(p):
    return list(csv.DictReader(p.open(),delimiter='\t'))
def tsv(rows,header=True):
    stream=io.StringIO()
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    if header:writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()
def queue():
    raw=subprocess.check_output(['squeue','-h','-r','-j','12812439','-o','%i|%T|%N|%r'],text=True)
    return raw,{p[0]:p[1:] for p in (line.split('|') for line in raw.splitlines())}
rows=read(prep/'package/designs.tsv')
assert len(rows)==25
groups={g:[r for r in rows if r['destination']==g] for g in ['fennario','compsci','DCC']}
assert [len(groups[g]) for g in groups]==[5,11,9]
with (out/'input_verification.log').open('w') as f:
    subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=prep/'package',stdout=f,stderr=subprocess.STDOUT,check=True)
old_dcc=read(campaign/'dcc_designs.tsv');old_local=read(campaign/'local_designs.tsv')
assert len(old_dcc)==91 and len(old_local)==48
carry=next(r for r in old_local if r['design_id']=='4wwi_flex_p1')
new_dcc=old_dcc+[carry]
new_local=[r for r in old_local if r['design_id']!=carry['design_id']]
assert len(new_dcc)==92 and len(new_local)==47
all139={r['design_id'] for r in new_dcc+new_local}
assert len(all139)==139
added24={r['design_id'] for r in rows if r['kind']!='existing_pending_transfer'}
assert len(added24)==24 and not added24&all139
raw,before=queue();(out/'queue_before.txt').write_text(raw)
carry_job='12812439_24'
assert before[carry_job][0]=='PENDING' and before[carry_job][2]=='Resources',before.get(carry_job)
assert not (Path('/usr/xtmp/lz280/markstar_local57_cpu64_20261005/A12812439')/carry['design_id']).exists()
old_files={name:(campaign/name).read_bytes() for name in ['dcc_cases.txt','dcc_designs.tsv','local_cases.txt','local_designs.tsv','summary.json']}
for name,data in old_files.items():(out/('before_'+name)).write_bytes(data)
submitted={};held=False;released=False
try:
    subprocess.run(['scontrol','hold',carry_job],check=True)
    held=True
    raw,after=queue();(out/'queue_after_hold.txt').write_text(raw)
    assert after[carry_job][0]=='PENDING' and after[carry_job][2]=='JobHeldUser'
    for group in ['fennario','compsci']:
        count=len(groups[group]);concurrency=2 if group=='fennario' else 4
        args=['sbatch','--parsable','--hold','--account=grisman',
              '--job-name=markstar-ext-'+group,'--partition='+('grisman' if group=='fennario' else 'compsci'),
              '--array=0-'+str(count-1)+'%'+str(concurrency),
              '--export=ALL,MARKSTAR_EXTENSION_PREP='+str(prep)+',MARKSTAR_EXTENSION_GROUP='+group]
        if group=='fennario':args+=['--constraint=a5000','--exclude=grisman-37,grisman-40,jerry[1-7]']
        else:args+=['--exclude=linux[31-40]']
        args.append(str(prep/'markstar_array.slurm'))
        job=subprocess.check_output(args,text=True).strip().split(';')[0]
        assert job.isdigit(),job
        submitted[group]=dict(array_job=job,count=count,concurrency=concurrency,cpus=64,memory_gib=192,limit_hours=336)
        (out/'submitted.json').write_text(json.dumps(submitted,indent=2)+'\n')
        details=subprocess.check_output(['scontrol','show','job',job,'-o'],text=True)
        (out/(group+'_held_job.txt')).write_text(details)
        assert 'Account=grisman' in details and 'CPUs/Task=64' in details and 'TimeLimit=14-00:00:00' in details
        assert 'JobHeldUser' in details
    # Append the additional completed-history case without renumbering any prior DCC row.
    (campaign/'dcc_cases.txt').write_bytes(old_files['dcc_cases.txt']+(carry['design_id']+'\n').encode())
    (campaign/'dcc_designs.tsv').write_bytes(old_files['dcc_designs.tsv']+tsv([carry],False).encode())
    (campaign/'local_cases.txt').write_text(''.join(r['design_id']+'\n' for r in new_local))
    (campaign/'local_designs.tsv').write_text(tsv(new_local))
    old_summary=json.loads(old_files['summary.json'])
    summary=dict(old_summary)
    summary.update(dcc_count=92,local_count=47,
                   total_hours=sum(float(r['baseline_mark_seconds']) for r in new_dcc)/3600,
                   longest_hours=max(float(r['baseline_mark_seconds']) for r in new_dcc)/3600,
                   local_total_hours=sum(float(r['baseline_mark_seconds']) for r in new_local)/3600,
                   extension_handoff_completed_case='4wwi_flex_p1',
                   previous_counts=dict(dcc=91,local=48),
                   basis='Historical timing provenance; the original 91 DCC rows retain their order. One unstarted 4wwi tier appended to keep its extension on one CPU model.')
    (campaign/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    records=[]
    for r in rows:
        r=dict(r)
        if r['destination'] in submitted:
            job=submitted[r['destination']]['array_job']
            r.update(status='SUBMITTED_LOCAL',slurm_job_id=job+'_'+r['group_index'])
        else:
            r.update(status='AWAITING_DCC_SUBMISSION',slurm_job_id='NOT_SUBMITTED')
        records.append(r)
    write_local=[r for r in records if r['destination']!='DCC']
    handoff=[r for r in records if r['design_id']=='4wwi_flex_p1']+[r for r in records if r['destination']=='DCC' and r['design_id']!='4wwi_flex_p1']
    for i,r in enumerate(handoff):r['group_index']=str(i)
    (here/'markstar_local16_submitted.tsv').write_text(tsv(write_local))
    (here/'dcc_submit9_designs.tsv').write_text(tsv(handoff))
    (here/'dcc_submit9_cases.txt').write_text(''.join(r['design_id']+'\n' for r in handoff))
    (here/'markstar24_status.tsv').write_text(tsv([r for r in records if r['kind']!='existing_pending_transfer']))
    handoff_record=dict(date='2026-10-06',status='AWAITING_DCC_SUBMISSION',
                       user_instruction='Push the manifest for remote submission as before',
                       designs=[r['design_id'] for r in handoff],historical_noncompletion=8,existing_completed_history_transfer=1,
                       cpus=64,time_limit_hours=336,memory_gib=192,combined_dcc_concurrency=6,
                       initial_added_long_concurrency=3,held_local_copy=dict(design='4wwi_flex_p1',job=carry_job,status='PENDING_HELD'),
                       previously_transferred_nine='Already externally submitted; do not repeat',
                       retire_local_copy_only_after_external_acceptance=True,
                       cpu_model_rule='Keep all newly launched tiers of each system on one actual CPU model',
                       frozen_input_package=str(prep/'package'),launch_audit=str(out))
    (here/'dcc_handoff.json').write_text(json.dumps(handoff_record,indent=2)+'\n')
    registry.update(markstar_status='LOCAL_SUBMITTED_DCC_HANDOFF_PENDING',markstar_prep=str(prep),
                    markstar_arrays=submitted,markstar_local_count=16,markstar_dcc_waiting=8,
                    dcc_extra_transfer='4wwi_flex_p1',markstar_packstar_dependency=None,
                    markstar_launch_audit=str(out),markstar_active_registry=str(here/'markstar24_status.tsv'))
    (here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
    assert len(read(here/'markstar_local16_submitted.tsv'))==16
    assert len(read(here/'dcc_submit9_designs.tsv'))==9
    assert (campaign/'dcc_designs.tsv').read_bytes().startswith(old_files['dcc_designs.tsv'])
    assert len(all139|added24)==163
    subprocess.run(['git','diff','--check'],cwd=repo,check=True)
    subprocess.run(['scontrol','release',submitted['fennario']['array_job'],submitted['compsci']['array_job']],check=True)
    released=True
except BaseException:
    if not released:
        for item in submitted.values():
            subprocess.run(['scancel',item['array_job']],check=False)
        for name,data in old_files.items():(campaign/name).write_bytes(data)
        if held:subprocess.run(['scontrol','release',carry_job],check=False)
        registry['markstar_status']='SUBMISSION_ROLLED_BACK_REVIEW_AUDIT'
        (here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
    raise
result=dict(submitted=submitted,local_submitted=16,dcc_waiting=9,
            old139_counts=dict(DCC=92,local=47),all163_counts=dict(DCC=100,local=63),
            held_local_copy=carry_job,packstar_dependency=None,output=str(out))
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print('LAUNCHED',json.dumps(result),flush=True)
