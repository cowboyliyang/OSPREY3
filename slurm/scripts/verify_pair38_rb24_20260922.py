"""Read-only startup verification; run within Slurm after launcher completes."""
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

assert os.environ.get('SLURM_JOB_ID')
root=Path(os.environ['PAIR38_VERIFY_ROOT'])
dest=root/('startup_audit_'+os.environ['SLURM_JOB_ID'])
dest.mkdir()
submission=json.loads((root/'submission.json').read_text())
jobs=[json.loads(line) for line in (root/'jobs.jsonl').read_text().splitlines()]
assert len(jobs)==len({j['job'] for j in jobs})==77
assert len(submission['jobs']['2'])==len(submission['jobs']['4'])==38
assert Counter((j['rb'],j['arm']) for j in jobs)=={(2,'pair-only'):38,(4,'pair-only'):38,(4,'budget-forward'):1}
records=[]
errors=[]
total=0
for group_name in ('rb2','rb4','rb4_triple_4z80'):
    group=root/group_name
    plan=json.loads((group/'plan.json').read_text())
    for name,digest in plan['source_sha256'].items():
        assert hashlib.sha256((group/'source'/name).read_bytes()).hexdigest()==digest
    for case in plan['cases']:
        total+=case['expected_rows']
        assert len(case['expected_sequences'])==len(set(case['expected_sequences']))==case['expected_rows']
        assert 0<case['host_gib']<case['heap_gib']<case['mem_gib']
        job=next(j for j in jobs if j['rb']==plan['rb'] and j['arm']==plan['arm'] and j['design']==case['design'])
        assert '--account=grisman' in job['command']
        matches=list((group/'runs').glob(case['design']+'_*/manifest.json'))
        assert len(matches)<=1
        record=dict(rb=plan['rb'],arm=plan['arm'],design=case['design'],job=job['job'],manifest_exists=bool(matches))
        if matches:
            path=matches[0]
            m=json.loads(path.read_text())
            p=m['properties']
            assert p['branchdp.cutoff.residualBudget']==str(plan['rb'])
            assert p['packstar.pac.frequencySeverity.tripleEta']==('false' if plan['arm']=='pair-only' else 'true')
            assert p['packstar.pac.frequencySeverity.jointMomentLearning']=='true'
            assert p['packstar.pac.frequencySeverity.proposalLearning']=='true'
            assert m['options']['seed']==42 and m['options']['arm']==plan['arm']
            assert m['pdb_sha256']==case['pdb_sha256']
            assert m['cpus']==case['cpus']
            record.update(status=m['status'],node=m['node'])
            if m.get('exit_code',0)!=0:
                errors.append(dict(record,error='nonzero Java exit',exit_code=m['exit_code']))
            log=path.parent/'run.log'
            if log.exists():
                with log.open(errors='replace') as stream:
                    start=stream.read(262144)
                for pattern in ('Exception in thread','OutOfMemoryError','CUDA_ERROR','UnsatisfiedLinkError'):
                    if pattern in start:
                        errors.append(dict(record,error=pattern))
        task_path=group/'tasks'/(str(case['index'])+'.json')
        if task_path.exists():
            task=json.loads(task_path.read_text())
            if task.get('status') in ('PROCESS_FAILED','AUDIT_FAILED'):
                errors.append(dict(record,error=task.get('error',task['status'])))
        records.append(record)
assert total==2245
scheduler=subprocess.check_output(['squeue','-h','-j',','.join(j['job'] for j in jobs),'-o','%i|%T|%P|%C|%m|%R'],text=True)
(dest/'squeue.txt').write_text(scheduler)
accounting=subprocess.check_output(['sacct','-X','-n','-P','-j',','.join(j['job'] for j in jobs),
    '--format=JobIDRaw,State,ExitCode,Elapsed,Partition,Account,AllocCPUS,ReqMem'],text=True)
(dest/'sacct.txt').write_text(accounting)
states=Counter()
for line in accounting.splitlines():
    f=line.split('|')
    assert f[5]=='grisman'
    states[f[1]]+=1
    if f[1] in ('FAILED','OUT_OF_MEMORY','NODE_FAIL','CANCELLED','TIMEOUT'):
        errors.append(dict(job=f[0],error='scheduler '+f[1]))
result=dict(root=str(root),planned_tasks=77,planned_sequences=total,manifests_checked=sum(r['manifest_exists'] for r in records),
    states=dict(states),errors=errors,groups={g:sum(r['manifest_exists'] for r in records if r['rb']==rb and r['arm']==arm)
        for g,rb,arm in [('rb2',2,'pair-only'),('rb4',4,'pair-only'),('rb4_triple_4z80',4,'budget-forward')]})
(dest/'records.json').write_text(json.dumps(records,indent=2)+'\n')
(dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
assert not errors, errors
