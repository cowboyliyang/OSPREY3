"""Freeze exact MARK* inputs for 16 local jobs and a nine-design DCC handoff."""
import csv
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
root=Path('/usr/xtmp/lz280/frontier_extension7_20261006')/('markstar_prep_'+os.environ['SLURM_JOB_ID'])
root.mkdir(parents=True,exist_ok=False)
package=root/'package'
base=Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
new_package=Path('/usr/xtmp/lz280/frontier_extension7_20261006/prep_12826503/package')
def read(p):
    return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
for name in ['markstar_runner.py','launch_markstar.py']:
    compile((here/name).read_text(),str(here/name),'exec')
subprocess.run(['bash','-n',str(here/'markstar_array.slurm')],check=True)
plan=read(here/'markstar24_proposed.tsv')
assert len(plan)==24 and len({r['design_id'] for r in plan})==24
specs={r['design_id']:r for r in read(base/'designs.tsv')}
new={r['design_id']:r for r in read(here/'new_designs.tsv')}
paper=json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']
runs={r['design_id']:r for r in paper['runs'] if r['method']=='MARK*'}
completed={r['design_id']:r for r in paper['pairs']}
bounds={r['design_id']:r for r in paper['censored_pairs']}
assert {r['design_id'] for r in plan if r['kind']=='historical_14day_noncompletion'}==set(bounds)
plan.append(dict(design_id='4wwi_flex_p1',system='4wwi',kind='existing_pending_transfer',destination='DCC',requested_cpus='64',limit_hours='336',state='PENDING_LOCAL_COPY_TO_HOLD'))
sources={r['system']:(new_package if r['design_id'] in new else base)/(new.get(r['design_id']) or specs[r['design_id']])['pdb_relative'] for r in plan}
print('STORAGE_ESTIMATE',json.dumps(dict(pdb_files=len(sources),pdb_bytes=sum(p.stat().st_size for p in sources.values()),
                                      preparation_files_under=250,preparation_bytes_under=100*2**20,
                                      run_outputs='Sixteen fresh run directories, caches, and partial scores under xtmp; no previous output copied.')),flush=True)
for name in ['structures','preflight']:
    (package/name).mkdir(parents=True)
for p in sources.values():
    shutil.copy2(p,package/'structures'/p.name)
for name in ['markstar_runner.py','markstar_array.slurm','launch_markstar.py']:
    shutil.copy2(here/name,root/name)
shutil.copy2(__file__,root/'prepare_markstar.py')
rows=[]; configs={}; group_index={'fennario':0,'compsci':0,'DCC':0}
for item in plan:
    name=item['design_id']; spec=new.get(name) or specs[name]
    reference=spec['anchor_design'] if name in new else name
    old=Path(runs[reference]['path'])
    command=shlex.split((old/'command.sh').read_text())
    assert '-Dosprey.bench.method=markstar' in command
    assert '-Dosprey.bench.numCPUs=16' in command
    cp=command[command.index('-cp')+1]
    assert all(Path(p).exists() for p in cp.split(':'))
    replacements={'-XX:ActiveProcessorCount=':'64','-Dosprey.bench.numCPUs=':'64',
                  '-Dosprey.bench.designId=':name,'-Dosprey.bench.outputDir=':'{RUN}',
                  '-Djava.io.tmpdir=':'{RUN}/tmp','-Dosprey.bench.pdbPath=':str(package/spec['pdb_relative']),
                  '-Dosprey.bench.mutable=':spec['mutable'],'-Dosprey.bench.flexible=':spec['flexible']}
    seen=set(); rewritten=[]
    for token in command:
        matched=next((key for key in replacements if token.startswith(key)),None)
        if matched:
            seen.add(matched);token=matched+replacements[matched]
        rewritten.append(token.replace(str(old),'{RUN}'))
    assert seen==set(replacements),(name,seen)
    heap=next(t for t in rewritten if t.startswith('-Xmx'))
    match=re.fullmatch(r'-Xmx(\d+)g',heap)
    assert match and int(match[1])<192,heap
    assert sha(package/spec['pdb_relative'])==spec['pdb_sha256']
    group=item['destination'];index=group_index[group];group_index[group]+=1
    row=dict(design_id=name,system=item['system'],kind=item['kind'],destination=group,group_index=index,
             flex_delta=spec['flex_delta'],mutable=spec['mutable'],flexible=spec['flexible'],
             total_positions=spec['total_positions'],expected_sequences=spec['expected_sequences'],
             pdb_relative=spec['pdb_relative'],pdb_sha256=spec['pdb_sha256'],requested_cpus=64,
             requested_memory_gib=192,java_heap_gib=int(match[1]),limit_hours=336,
             historical_mark_seconds=completed[name]['mark_seconds'] if name in completed else '',
             historical_mark_lower_seconds=bounds[name]['mark_lower_seconds'] if name in bounds else '',
             reference_design=reference,reference_mark_path=str(old),
             source_package=str(new_package if name in new else base),
             status='PREPARED_NOT_SUBMITTED')
    rows.append(row)
    configs[name]=dict(row=row,command_template=rewritten,source_command_sha256=sha(old/'command.sh'),
                       inherited_protocol=json.loads((old/'protocol.json').read_text()) if (old/'protocol.json').exists() else None,
                       model_contains='Xeon(R) Gold 5320' if group=='fennario' else 'AMD EPYC 9554' if group=='compsci' else None,
                       input_pdb=str(package/spec['pdb_relative']))
    preflight=package/'preflight'/name;preflight.mkdir()
    precommand=[]
    for token in rewritten:
        token=token.replace('{RUN}',str(preflight))
        if token.startswith('-Xmx'):token='-Xmx8g'
        elif token.startswith('-Xms'):token='-Xms1g'
        elif token.startswith('-XX:ActiveProcessorCount='):token='-XX:ActiveProcessorCount=4'
        elif token.startswith('-Dosprey.bench.numCPUs='):token='-Dosprey.bench.numCPUs=4'
        elif token.startswith('-Dosprey.bench.method='):token='-Dosprey.bench.method=sequence_dump'
        precommand.append(token)
    (preflight/'tmp').mkdir()
    with (preflight/'preflight.log').open('w') as f:
        subprocess.run(precommand,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=240)
    log=(preflight/'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log,name
    positions=re.findall(r'\[FRONTIER_POSITION\] state=Complex residue=(\S+)',log)
    assert len(positions)==int(spec['total_positions'])
    assert set(positions)==set((spec['mutable']+';'+spec['flexible']).split(';'))
    sequence_file=name+'_sequences.tsv'
    expected=read((new_package if name in new else base)/'preflight'/name/sequence_file)
    actual=read(preflight/sequence_file)
    assert actual==expected and len(actual)==int(spec['expected_sequences']),name
    print('PREFLIGHT_OK',name,'group',group,'positions',len(positions),'sequences',len(actual),'heap',match[1],flush=True)
assert group_index=={'fennario':5,'compsci':11,'DCC':9}
write(package/'designs.tsv',rows)
(package/'run_configs.json').write_text(json.dumps(configs,indent=2)+'\n')
(package/'protocol.json').write_text(json.dumps(dict(cpus=64,time_limit_hours=336,memory_gib=192,
    heap='Preserved from each historical reference command',epsilon=0.683,stability_filter=False,
    hardware={'fennario':'Xeon Gold 5320','compsci':'AMD EPYC 9554'},
    packstar_dependency=None,user_authorization='Run now; DCC submission through pushed manifests as before',
    counts=group_index),indent=2)+'\n')
checks=[sha(p)+'  '+str(p.relative_to(package)) for p in sorted(package.rglob('*')) if p.is_file()]
(package/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(root/'verification.json').write_text(json.dumps(dict(preflight_passed=25,local=16,dcc_handoff=9,
    original_completed_transfer=1,new_plus_timeout=24,counts=group_index),indent=2)+'\n')
(root/'READY').write_text('25 actual position and sequence preflights passed\n')
print('READY',root,flush=True)
