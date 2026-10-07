"""Retry failed GPU groups with explicit memory, preserving frozen inputs."""
import json,os,shutil,subprocess
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
here=Path(__file__).parent
registry=here/'registry.json'
assert not registry.exists(), 'Already submitted; inspect registry before any retry'
root=Path('/usr/xtmp/lz280/packstar_gpu_memory_retry_20261006')/('prep_'+os.environ['SLURM_JOB_ID'])
root.mkdir(parents=True,exist_ok=False)
specs=[('deadline_focus_20261006','/usr/xtmp/lz280/markstar_deadline_focus_20261006/pack_pairs_12827168','12827169'),
       ('near14_20261006','/usr/xtmp/lz280/markstar_near14_20261006/prep_12827318','12827321')]
record=dict(reason='Explicit 1000 GiB request prevents inherited small SLURM_MEM_PER_NODE with mem=0',
            account='grisman',cpus=128,gpus=4,memory_gib=1000,limit_hours_per_group=6,
            markstar_changes=[],groups=[],state='PREPARING')
def save():registry.write_text(json.dumps(record,indent=2)+'\n')
save()
previous=None
for campaign,oldpath,oldjob in specs:
    state=subprocess.check_output(['sacct','-n','-X','-j',oldjob,'--format=State','-P'],text=True).strip()
    assert state.startswith('FAILED'),(oldjob,state)
    old=Path(oldpath);new=root/campaign;new.mkdir()
    for name in ['package','build']:(new/name).symlink_to(old/name,target_is_directory=True)
    for name in ['READY','portable_run.py','verification.json','run_group.py']:
        shutil.copy2(old/name,new/name)
    shutil.copy2(here.parent/campaign/'run_group.slurm',new/'run_group.slurm')
    subprocess.run(['bash','-n',str(new/'run_group.slurm')],check=True)
    if campaign.startswith('deadline'):
        os.environ['FOCUS_PACK_INDICES']=','.join(str(i) for i in range(12))
    args=['sbatch','--parsable','--hold','--account=grisman','--mem=1000G',
          '--export=ALL,EXTENSION_PREP='+str(new)+',SLURM_MEM_PER_NODE=1024000']
    if previous:args.append('--dependency=afterany:'+previous)
    job=subprocess.check_output([*args,str(new/'run_group.slurm')],text=True).strip().split(';')[0]
    assert job.isdigit()
    output=Path('/usr/xtmp/lz280')/('markstar_deadline_focus_20261006' if campaign.startswith('deadline') else 'markstar_near14_20261006')/('pack_'+job)
    group=dict(campaign=campaign,old_job=oldjob,job=job,prep=str(new),output=str(output),after_job=previous)
    record['groups'].append(group);save()
    detail=subprocess.check_output(['scontrol','show','job',job,'-o'],text=True)
    assert 'MinMemoryNode=1000G ' in detail or 'mem=1000G' in detail or 'MinMemoryNode=1024000M ' in detail,detail
    (new/'submitted_job.txt').write_text(detail)
    phase=here.parent/campaign
    if campaign.startswith('deadline'):
        for name in ['screen_pack.py','screen_pack.slurm']:shutil.copy2(old/name,new/name)
        auditargs=['--export=ALL,EXTENSION_PREP='+str(new)+',FOCUS_PACK_ROOT='+str(output)+',FOCUS_PACK_JOB='+job,str(new/'screen_pack.slurm')]
    else:
        auditargs=['--export=ALL,PACK_RUN_ROOT='+str(output),str(phase/'collect_pack.slurm')]
    audit=subprocess.check_output(['sbatch','--parsable','--account=grisman','--dependency=afterany:'+job,*auditargs],text=True).strip().split(';')[0]
    group['audit_job']=audit;save()
    active=json.loads((phase/'registry.json').read_text())
    active['gpu_memory_retry_registry']=str(registry)
    active['previous_failed_packstar_job']=oldjob
    active['packstar_job']=job;active['packstar_output']=str(output)
    active['packstar_after_job']=previous
    active['packstar_memory_gib']=1000
    active['packstar_retry_audit_job']=audit
    (phase/'registry.json').write_text(json.dumps(active,indent=2)+'\n')
    previous=job
for group in record['groups']:subprocess.run(['scontrol','release',group['job']],check=True)
record['state']='SUBMITTED';save()
print(json.dumps(record,indent=2),flush=True)
