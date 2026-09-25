"""Slurm-only validation, input audit, and six disjoint single-GPU arrays."""
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

job = os.environ['SLURM_JOB_ID']
repo = Path('/home/users/lz280/OSPREY3-fresh-packstar')
base = Path('/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056')
root = Path('/usr/xtmp/lz280/packstar_rb38_20260918') / f'launch_{job}'
root.mkdir(parents=True, exist_ok=False)
source = root / 'source'
source.mkdir()
assert (base / 'READY').is_file()
for name in ['run.py', 'baseline38.csv', 'production.properties']:
    origin = repo / 'slurm/h200' / name if name == 'run.py' else base / 'source/slurm/h200' / name
    shutil.copy2(origin, source / name)
shutil.copy2(repo / 'slurm/scripts/run_rb38_20260918.slurm', source / 'run.slurm')
shutil.copy2(__file__, source / 'launch.py')
ast.parse((source / 'run.py').read_text())
subprocess.run(['bash', '-n', str(source / 'run.slurm')], check=True)
with (source / 'baseline38.csv').open() as stream:
    rows = list(csv.reader(line for line in stream if not line.startswith('#')))
assert len(rows) == 38 and len({r[0] for r in rows}) == 38
inputs = []
for index, row in enumerate(rows):
    pdb = Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped') / row[1] / f'{row[1]}.min.reduce.renum.pdb'
    data = pdb.read_bytes()
    inputs.append(dict(index=index, system=row[0], bytes=len(data),
                       pdb=str(pdb), sha256=hashlib.sha256(data).hexdigest()))
audit = dict(build=str(base), inputs=inputs, rb=[1, 2, 4], seed=42,
             epsilon=0.683, tasks=114, gpu_count=1, gpu_type='unrestricted',
             cpus=8, memory_gib=112, heap_gib=96, host_budget_gib=80,
             gpu_budget='floor(80% of allocated VRAM in GiB)',
             storage_estimate='Source/input audit under 10 MiB; allow 100 GiB reports/EMAT and up to 57 TiB transient DP files across 114 tasks; approximately 100000 files, workload dependent.',
             source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()})
(root / 'manifest.json').write_text(json.dumps(audit, indent=2))
(root / 'READY').write_text('Inputs and submission scripts validated; existing verified Java/CUDA build reused.\n')
with (root / 'jobs.jsonl').open('x') as output:
    for rb in [1, 2, 4]:
        for side, partition, array in [('fennario', 'grisman', '0-36:2'), ('compsci-gpu', 'compsci-gpu', '1-37:2')]:
            command = ['sbatch', '--parsable', '--account=grisman', f'--partition={partition}',
                       f'--array={array}', f'--job-name=rb{rb}_38_{side}',
                       f'--export=ALL,SWEEP_ROOT={root},RB={rb}']
            if side == 'fennario':
                command += ['--nodelist=fennario-[01-06]']
            command += [str(source / 'run.slurm')]
            ident = subprocess.check_output(command, text=True).strip().split(';')[0]
            record = dict(rb=rb, side=side, array=array, job=ident, tasks=19)
            output.write(json.dumps(record) + '\n')
            output.flush()
            print(json.dumps(record), flush=True)
(root / 'SUBMITTED').write_text('114 tasks submitted in six arrays.\n')
print(root, flush=True)
