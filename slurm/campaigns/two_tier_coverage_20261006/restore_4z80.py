"""Restore only 4z80 p0 on matching hardware, before a longer continuation."""
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here=Path(__file__).parent
active=here.parent/'deadline_focus_20261006'
assert not (here/'registry.json').exists(), 'Already applied or partially submitted; inspect registry'
out=Path('/usr/xtmp/lz280/markstar_two_tier_20261006')/('prep_'+os.environ['SLURM_JOB_ID'])
package=out/'package'
package.mkdir(parents=True,exist_ok=False)
base=Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
pack=Path('/usr/xtmp/lz280/packstar_gpu101_20261005/run_12814879/results/4z80_flex_p0_full_pair-only_s42_J12814879_T143')
name='4z80_flex_p0'
parent='12812439_41'
child='12827032_2'

def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def write(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def detail(job):return subprocess.check_output(['scontrol','show','job',job,'-o'],text=True)
def pending(job):
    d=detail(job)
    assert 'JobState=PENDING ' in d and 'Restarts=0 ' in d and 'RunTime=00:00:00 ' in d,d
    return d

row=next(r for r in read(base/'designs.tsv') if r['design_id']==name)
pdb=base/row['pdb_relative']
print('STORAGE_ESTIMATE',json.dumps(dict(existing_structure_bytes=pdb.stat().st_size,
    preparation_files_under=30,preparation_bytes_under=2*2**20,
    production='Fresh MARK* matrices/logs under scratch; no dataset or build copy',root=str(out))),flush=True)
assert digest(pdb)==row['pdb_sha256']
pm=json.loads((pack/'manifest.json').read_text())
assert pm['exit_code']==0 and pm['result_rows']==pm['expected_rows']==39
assert pm['options']['gpus']==4
for key in ['mutable','flexible','expected_sequences','pdb_sha256']:
    assert row[key]==pm['design'][key],key
pack_rows=list(csv.DictReader((pack/(name+'_packstar.csv')).open()))
estimated=sum(all(r[k]=='Estimated' for k in ['prot_status','lig_status','comp_status']) for r in pack_rows)
assert len(pack_rows)==39 and estimated==33
pack_seconds=float(re.search(r'elapsed=([\d.]+)',(pack/'wall.time').read_text())[1])
anchor=Path('/usr/xtmp/lz280/markstar_local57_cpu64_20261005/A12812439/4z80_flex_p1')
assert 'AMD EPYC 9554' in (anchor/'lscpu.txt').read_text()
assert json.loads((anchor/'run_manifest.json').read_text())['cpus']==64
history=json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']['runs']
source=Path(next(r['path'] for r in history if r['method']=='MARK*' and r['design_id']==name))
command=shlex.split((source/'command.sh').read_text())
changes={'-XX:ActiveProcessorCount=':'64','-Dosprey.bench.numCPUs=':'64',
         '-Dosprey.bench.designId=':name,'-Dosprey.bench.outputDir=':'{RUN}',
         '-Djava.io.tmpdir=':'{RUN}/tmp','-Dosprey.bench.pdbPath=':str(pdb),
         '-Dosprey.bench.mutable=':row['mutable'],'-Dosprey.bench.flexible=':row['flexible']}
changed=set();rewritten=[]
for token in command:
    key=next((k for k in changes if token.startswith(k)),None)
    if key:changed.add(key);token=key+changes[key]
    rewritten.append(token.replace(str(source),'{RUN}'))
assert changed==set(changes) and '-Dosprey.bench.method=markstar' in rewritten
assert all(Path(p).exists() for p in rewritten[rewritten.index('-cp')+1].split(':'))
heap=int(re.fullmatch(r'-Xmx(\d+)g',next(t for t in rewritten if t.startswith('-Xmx')))[1])
assert heap<192
preflight=out/'preflight';preflight.mkdir();(preflight/'tmp').mkdir()
precommand=[]
for token in rewritten:
    token=token.replace('{RUN}',str(preflight))
    if token.startswith('-Xmx'):token='-Xmx8g'
    elif token.startswith('-XX:ActiveProcessorCount='):token='-XX:ActiveProcessorCount=4'
    elif token.startswith('-Dosprey.bench.numCPUs='):token='-Dosprey.bench.numCPUs=4'
    elif token.startswith('-Dosprey.bench.method='):token='-Dosprey.bench.method=sequence_dump'
    precommand.append(token)
with (preflight/'run.log').open('w') as f:
    subprocess.run(precommand,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=240)
log=(preflight/'run.log').read_text()
assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log
seq=read(preflight/(name+'_sequences.tsv'))
original=read(base/'preflight'/name/(name+'_sequences.tsv'))
assert len(seq)==39 and [r['sequence'] for r in seq]==[r['sequence'] for r in original]
assert {r['sequence'] for r in seq}=={r['sequence'] for r in pack_rows}
row.update(destination='compsci',group_index=0,requested_cpus=64,requested_memory_gib=192,
           java_heap_gib=heap,limit_hours=72,kind='second_tier_coverage',
           reference_design=name,reference_mark_path=str(source))
config=dict(row=row,command_template=rewritten,input_pdb=str(pdb),model_contains='AMD EPYC 9554',
            source_command_sha256=digest(source/'command.sh'))
write(package/'designs.tsv',[row])
(package/'run_configs.json').write_text(json.dumps({name:config},indent=2)+'\n')
for filename in ['markstar_runner.py','markstar_array.slurm']:
    if filename.endswith('.py'):compile((here/filename).read_text(),filename,'exec')
    else:subprocess.run(['bash','-n',str(here/filename)],check=True)
    shutil.copy2(here/filename,out/filename)
(out/'READY').write_text('MARK* preflight and existing GPU PACK* match all 39 sequences\n')
write(here/'designs.tsv',[row])

actions=read(active/'design_actions.tsv')
action=next(r for r in actions if r['design']==name)
assert action['action']=='DEFER_UNSTARTED' and action['current_status']=='LOCAL_PENDING_COPY_CANCELLED'
edges=read(active/'local_lane_dependencies.tsv')
edge=next(r for r in edges if r['job']==child)
assert edge['parent']==parent
assert 'afterany:'+parent in pending(child)
pending(parent)
(out/'child_before.txt').write_text(detail(child))
before_p1=detail('12812439_7');(out/'p1_before.txt').write_text(before_p1)
assert 'JobState=RUNNING ' in before_p1
args=['sbatch','--parsable','--account=grisman','--hold','--job-name=markstar-4z80-p0',
      '--partition=compsci','--exclude=linux[31-40]','--array=0-0%1',
      '--dependency=afterany:'+parent,'--export=ALL,MARKSTAR_EXTENSION_PREP='+str(out)+',MARKSTAR_EXTENSION_GROUP=compsci',
      str(out/'markstar_array.slurm')]
job=subprocess.check_output(args,text=True).strip().split(';')[0]
assert job.isdigit()
registry=dict(prep_job=os.environ['SLURM_JOB_ID'],prep=str(out),design=name,job=job+'_0',
    after_job=parent,before_job=child,cpus=64,memory_gib=192,limit_hours=72,cpu_model='AMD EPYC 9554',
    packstar_existing_output=str(pack),packstar_seconds=pack_seconds,packstar_rows=39,packstar_fully_estimated=estimated,
    new_packstar_jobs=[],preflight='PASS',state='SUBMITTED_HELD_AWAITING_CHAIN_INSERTION',dcc_mutations=[],
    user_decision='Only p0 and the existing running p1; no higher exploration or final-round target for 4z80')
def save(): (here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
save()
pending(child)
subprocess.run(['scontrol','hold',child],check=True)
pending(child)
subprocess.run(['scontrol','update','JobId='+child,'Dependency=afterany:'+job+'_0'],check=True)
assert 'afterany:'+job+'_0' in detail(child)
subprocess.run(['scontrol','release',child],check=True)
subprocess.run(['scontrol','release',job],check=True)
edge.update(parent=job+'_0',action='SECOND_TIER_COVERAGE_BEFORE_LONG_CONTINUATION')
edges.append(dict(job=job+'_0',parent=parent,action='RESTORED_SECOND_TIER_COVERAGE'))
assert len({r['parent'] for r in edges})==len(edges), 'Do not branch a serial CPU lane'
write(active/'local_lane_dependencies.tsv',edges)
action.update(action='KEEP_SECOND_TIER_COVERAGE',proposed_destination='compsci',current_local_job=job+'_0',current_status='SUBMITTED_CPU_ONLY')
write(active/'design_actions.tsv',actions)
coverage=read(active/'system_coverage.tsv')
c=next(r for r in coverage if r['system']=='4z80')
c.update(pending_or_new=name,fallback_until_completed=name)
write(active/'system_coverage.tsv',coverage)
policy=json.loads((active/'policy.json').read_text())
policy['two_tier_coverage_goal']=dict(prefer_two_observed_tiers_per_system_before_final_round=True,
    queued_is_not_completed=True,minimum_one_completed_under_14days=True,planned_systems_with_at_least_two_tiers=38,
    four_z80=dict(keep=['4z80_flex_p0','4z80_flex_p1'],higher_exploration=False,final_round_candidate=False,registry=str(here/'registry.json')))
policy['rules'].append('User chose 4z80 p0 as the second tier; retain running p1 and stop further 4z80 expansion. Prefer two observed tiers per system before the final long-run round.')
(active/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
old_registry=json.loads((active/'registry.json').read_text())
old_registry['two_tier_coverage_registry']=str(here/'registry.json')
(active/'registry.json').write_text(json.dumps(old_registry,indent=2)+'\n')
counts={s:sum(r['system']==s and not r['action'].startswith(('DEFER','HOLD')) for r in actions) for s in {r['system'] for r in actions}}
assert len(counts)==38 and min(counts.values())>=2
submitted=detail(job+'_0')
assert 'TimeLimit=3-00:00:00 ' in submitted and 'CPUs/Task=64 ' in submitted and 'JobHeldUser' not in submitted
after_p1=detail('12812439_7')
assert re.search(r' StartTime=(\S+)',before_p1)[1]==re.search(r' StartTime=(\S+)',after_p1)[1]
assert re.search(r' TimeLimit=(\S+)',before_p1)[1]==re.search(r' TimeLimit=(\S+)',after_p1)[1]
registry.update(state='SUBMITTED',applied_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    all_38_systems_have_two_active_or_planned_tiers=True,original_running_task_unchanged=True,
    changed_existing_dependency=dict(job=child,old_parent=parent,new_parent=job+'_0'),
    command_sha256=digest(package/'run_configs.json'))
save();(out/'summary.json').write_text(json.dumps(registry,indent=2)+'\n')
subprocess.run(['git','diff','--check'],cwd=here.parents[2],check=True)
print('READY',json.dumps(registry),flush=True)
