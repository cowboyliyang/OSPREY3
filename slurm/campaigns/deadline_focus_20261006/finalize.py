"""Freeze MARK* commands, replace pending PACK* screening, and publish exact handoff."""
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
repo=here.parents[2]
registry=json.loads((here/'registry.json').read_text())
prep=Path(registry['prep_root'])
assert (prep/'READY').exists()
out=prep.parent/('launch_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rs):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def query():
    raw=subprocess.check_output(['squeue','-u','lz280','-r','-h','-o','%i|%T|%R'],text=True)
    return {p[0]:p[1:] for line in raw.splitlines() if (p:=line.split('|'))}
def submit(args):
    jid=subprocess.check_output(['sbatch','--parsable','--account=grisman',*args],text=True).strip().split(';')[0]
    assert jid.isdigit(),jid
    with (out/'submitted_jobs.txt').open('a') as f:f.write(jid+' '+shlex.join(args)+'\n')
    return jid
rows=read(prep/'package/designs.tsv')
new=[r for r in rows if r['kind']!='retained_extension']
assert len(new)==6
oldrows={r['design_id']:r for r in read(here.parent/'frontier_extension7_20261006/new_designs.tsv')}
for r in rows:
    if r['kind']=='retained_extension':
        for key in ['mutable','flexible','pdb_sha256','expected_sequences','total_positions']:
            assert r[key]==oldrows[r['design_id']][key],(r['design_id'],key)
paper=json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']
historical={r['design_id']:r for r in paper['runs'] if r['method']=='MARK*'}
markprep=prep/'markstar';package=markprep/'package';package.mkdir(parents=True)
configs={};markrows=[];indices={'fennario':0,'compsci':0,'DCC':0}
for r in new:
    source=Path(historical[r['anchor_design']]['path'])
    command=shlex.split((source/'command.sh').read_text())
    replace={'-XX:ActiveProcessorCount=':'64','-Dosprey.bench.numCPUs=':'64',
             '-Dosprey.bench.designId=':r['design_id'],'-Dosprey.bench.outputDir=':'{RUN}',
             '-Djava.io.tmpdir=':'{RUN}/tmp','-Dosprey.bench.pdbPath=':str(prep/'package'/r['pdb_relative']),
             '-Dosprey.bench.mutable=':r['mutable'],'-Dosprey.bench.flexible=':r['flexible']}
    changed=set();rewritten=[]
    for token in command:
        key=next((k for k in replace if token.startswith(k)),None)
        if key:changed.add(key);token=key+replace[key]
        rewritten.append(token.replace(str(source),'{RUN}'))
    assert changed==set(replace)
    assert '-Dosprey.bench.method=markstar' in rewritten
    assert all(Path(p).exists() for p in rewritten[rewritten.index('-cp')+1].split(':'))
    heap=int(re.fullmatch(r'-Xmx(\d+)g',next(t for t in rewritten if t.startswith('-Xmx')))[1])
    group=r['destination'];index=indices[group];indices[group]+=1
    row=dict(r,group_index=index,requested_cpus=64,requested_memory_gib=192,java_heap_gib=heap,limit_hours=336,
             reference_design=r['anchor_design'],reference_mark_path=str(source),source_package=str(prep/'package'))
    assert heap<192
    markrows.append(row)
    configs[r['design_id']]=dict(row=row,command_template=rewritten,source_command_sha256=hashlib.sha256((source/'command.sh').read_bytes()).hexdigest(),
        model_contains='Xeon(R) Gold 5320' if group=='fennario' else 'AMD EPYC 9554' if group=='compsci' else None,
        input_pdb=str(prep/'package'/r['pdb_relative']))
write(package/'designs.tsv',markrows)
(package/'run_configs.json').write_text(json.dumps(configs,indent=2)+'\n')
for filename in ['markstar_runner.py','markstar_array.slurm']:
    if filename.endswith('.py'):compile((here/filename).read_text(),filename,'exec')
    else:subprocess.run(['bash','-n',str(here/filename)],check=True)
    shutil.copy2(here/filename,markprep/filename)
(markprep/'READY').write_text('Six new definitions passed MARK* and PACK* position/sequence preflights\n')
for filename in ['screen_pack.py','screen_pack.slurm']:
    shutil.copy2(here/filename,prep/filename)
compile((here/'screen_pack.py').read_text(),'screen_pack.py','exec')
# Preserve every started GPU job. Only replace the old seven-case job while pending.
live=query();oldpack='12826506';replace_pack=oldpack in live and live[oldpack][0]=='PENDING'
if replace_pack:
    subprocess.run(['scontrol','hold',oldpack],check=True)
    live=query();replace_pack=oldpack in live and live[oldpack][0]=='PENDING'
selection=','.join(str(i) for i in range(10 if replace_pack else 6))
# Pass the comma-containing selection through the environment, not --export syntax.
os.environ['FOCUS_PACK_INDICES']=selection
pack=submit(['--hold','--export=ALL,EXTENSION_PREP='+str(prep),str(prep/'run_group.slurm')])
if replace_pack:
    subprocess.run(['scancel','--state=PENDING',oldpack],check=True)
    assert oldpack not in query() or query()[oldpack][0]!='PENDING'
subprocess.run(['scontrol','release',pack],check=True)
packroot=prep.parent/('pack_'+pack)
exports='ALL,EXTENSION_PREP='+str(prep)+',FOCUS_PACK_ROOT='+str(packroot)+',FOCUS_PACK_JOB='+pack
screen=submit(['--dependency=afterany:'+pack,'--export='+exports,str(prep/'screen_pack.slurm')])
# Eight compsci lanes and four fennario lanes, with coverage before one extra target.
dependencies={
 '12812439_52':'12812439_12','12812439_41':'12812439_7','12812439_26':'12812439_13',
 '12826674_0':'12812439_8','12826674_1':'12812439_4','12826674_3':'12812439_5',
 '12826674_7':'12812439_9','12826674_9':'12812439_11',
 '12826673_0':'12812439_1','12826673_3':'12812439_3','12826673_4':'12812439_0'}
lane_actions=[]
for child,parent in dependencies.items():
    live=query()
    if child not in live or live[child][0]!='PENDING':
        lane_actions.append(dict(job=child,parent=parent,action='SKIP_STARTED_OR_FINISHED'));continue
    # afterany tolerates a predecessor that ended without a usable result.
    subprocess.run(['scontrol','update','JobId='+child,'Dependency=afterany:'+parent],check=True)
    lane_actions.append(dict(job=child,parent=parent,action='PENDING_DEPENDENCY_SET'))
newjobs=[]
for group,parent in [('fennario','12812439_2'),('compsci','12812439_26')]:
    args=['--job-name=markstar-focus-'+group,'--partition='+('grisman' if group=='fennario' else 'compsci'),
          '--array=0-0%1','--dependency=afterok:'+screen+',afterany:'+parent,
          '--export=ALL,MARKSTAR_EXTENSION_PREP='+str(markprep)+',MARKSTAR_EXTENSION_GROUP='+group+',FOCUS_PACK_ROOT='+str(packroot)]
    if group=='fennario':args+=['--constraint=a5000','--exclude=grisman-37,grisman-40,jerry[1-7]']
    else:args+=['--exclude=linux[31-40]']
    jid=submit(args+[str(markprep/'markstar_array.slurm')])
    row=next(r for r in markrows if r['destination']==group)
    newjobs.append(dict(design=row['design_id'],job=jid+'_0',group=group,after_job=parent,screening_job=screen,
                        gate='PACK* normal completion <=1800 seconds; start by October 16'))
write(here/'new_local_submissions.tsv',newjobs)
write(here/'local_lane_dependencies.tsv',lane_actions)
write(here/'dcc_new_designs.tsv',[r for r in markrows if r['destination']=='DCC'])
(here/'dcc_new_cases.txt').write_text(''.join(r['design_id']+'\n' for r in markrows if r['destination']=='DCC'))
registry.update(status='LOCAL_ADJUSTED_REMOTE_HANDOFF_READY',packstar_status='SUBMITTED',markstar_status='TWO_NEW_LOCAL_JOBS_SUBMITTED_WITH_PACKSTAR_GATE',packstar_job=pack,packstar_replaced_pending_job=oldpack if replace_pack else None,
    packstar_selected_indices=selection,packstar_output=str(packroot),packstar_screening_job=screen,
    markstar_prep=str(markprep),new_local_jobs=newjobs,launch_audit=str(out),dcc_status='NOT_APPLIED_REMOTELY')
(here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
(out/'summary.json').write_text(json.dumps(registry,indent=2)+'\n')
print('READY',json.dumps(registry),flush=True)
