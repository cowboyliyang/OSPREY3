"""Run one preflighted DCC coverage design on 64 AMD EPYC compsci CPUs."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time

assert os.environ.get('SLURM_JOB_ID')
assert int(os.environ['SLURM_CPUS_PER_TASK'])==64
prep=Path(os.environ['MARKSTAR_EXTENSION_PREP'])
group=os.environ['MARKSTAR_EXTENSION_GROUP']
assert group in ('fennario','compsci') and (prep/'READY').is_file()
package=prep/'package'
rows=list(csv.DictReader((package/'designs.tsv').open(),delimiter='\t'))
index=int(os.environ['SLURM_ARRAY_TASK_ID'])
row=next(r for r in rows if r['destination']==group and int(r['group_index'])==index)
name=row['design_id']
config=json.loads((package/'run_configs.json').read_text())[name]
# MARK* and PACK* run independently. Numerical speedups require both results,
# but GPU queueing or a slow PACK* measurement must not block a CPU allocation.
node=os.environ['SLURMD_NODENAME']
if group=='fennario':
    assert re.fullmatch(r'fennario-0[1-6]',node),node
else:
    assert re.fullmatch(r'compsci-cluster-fitz-(3[5-9]|4[0-4])',node),node
cpu=subprocess.check_output(['lscpu'],text=True,env=dict(os.environ,LC_ALL='C'))
assert config['model_contains'] in cpu,(name,node,config['model_contains'])
assert hashlib.sha256(Path(config['input_pdb']).read_bytes()).hexdigest()==row['pdb_sha256']
root=Path('/usr/xtmp/lz280/markstar_dcc_coverage_to_compsci_20261009')/('markstar_A'+os.environ['SLURM_ARRAY_JOB_ID'])/name
root.mkdir(parents=True,exist_ok=False);(root/'tmp').mkdir()
command=[token.replace('{RUN}',str(root)) for token in config['command_template']]
assert '-Dosprey.bench.numCPUs=64' in command and '-XX:ActiveProcessorCount=64' in command
assert '-Dosprey.bench.method=markstar' in command
assert all(Path(p).exists() for p in command[command.index('-cp')+1].split(':'))
(root/'command.sh').write_text(shlex.join(command)+'\n')
(root/'lscpu.txt').write_text(cpu)
(root/'design.json').write_text(json.dumps(row,indent=2)+'\n')
(root/'protocol.json').write_text(json.dumps(config,indent=2)+'\n')
shutil.copy2(__file__,root/'runner.py')
with (root/'slurm_job.txt').open('w') as f:
    subprocess.run(['scontrol','show','job',os.environ['SLURM_JOB_ID']],stdout=f,stderr=subprocess.STDOUT,check=True)
manifest=dict(design_id=name,kind=row['kind'],group=group,node=node,cpus=64,
              memory_gib=192,heap_gib=int(row['java_heap_gib']),time_limit_hours=336,
              stability_filter=False,epsilon=0.683,precision='FP64',ccd='CPU',
              affinity=sorted(os.sched_getaffinity(0)),status='RUNNING',start_epoch=time.time(),
              slurm_job_id=os.environ['SLURM_JOB_ID'],array_job_id=os.environ['SLURM_ARRAY_JOB_ID'],
              array_task_id=index,prep=str(prep),reference_design=row['reference_design'],
              reference_mark_path=row['reference_mark_path'],source_command_sha256=config['source_command_sha256'])
(root/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('START',json.dumps(manifest),flush=True)
with (root/'run.log').open('w') as f:
    result=subprocess.run(['/usr/bin/time','-o',str(root/'wall.time'),'-f',
                           'elapsed=%e user=%U system=%S maxRssKiB=%M',*command],
                          cwd=root,stdout=f,stderr=subprocess.STDOUT)
output=root/(name+'_markstar.csv')
results=list(csv.DictReader(output.open())) if output.exists() else []
complete=result.returncode==0 and len(results)==int(row['expected_sequences'])
manifest.update(end_epoch=time.time(),exit_code=result.returncode,result_rows=len(results),
                expected_rows=int(row['expected_sequences']),status='COMPLETED' if complete else 'FAILED')
(root/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('END',json.dumps(manifest),flush=True)
raise SystemExit(0 if complete else result.returncode or 5)
