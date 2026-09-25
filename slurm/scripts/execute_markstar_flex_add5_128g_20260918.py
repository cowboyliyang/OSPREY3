"""Run the two pending MARK* designs with 128 GiB RAM and a 96 GiB heap."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

assert os.environ.get('SLURM_JOB_ID'), 'Run via Slurm'
prep, root = map(Path, sys.argv[1:3])
task = int(sys.argv[3])
assert task in (24, 25)
assert int(os.environ['SLURM_CPUS_PER_TASK']) == 16
assert int(os.environ['SLURM_MEM_PER_NODE']) == 128 * 1024
designs = json.loads((prep / 'package/designs.json').read_text())
row = designs[task]
assert row['task_id'] == task
assert row['design_id'] == f'5dc0_flex_p{task - 18}'
pdb = prep / 'package' / row['pdb_relative']
assert hashlib.sha256(pdb.read_bytes()).hexdigest() == row['pdb_sha256']
cp = (prep / 'test_classpath.txt').read_text().strip()
java = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'
cmd = [java,
       '--add-opens', 'java.base/java.util=ALL-UNNAMED',
       '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
       '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED',
       '-Xmx96g', '-Xms4g', '-XX:-UseSuperWord', '-XX:ActiveProcessorCount=16',
       '-Djava.io.tmpdir=' + str(root / 'tmp'), '-Dosprey.bench.method=markstar',
       '-Dosprey.bench.epsilon=0.683', '-Dosprey.bench.numCPUs=16', '-Dosprey.wmb.numGpus=0',
       '-Dosprey.bench.designId=' + row['design_id'], '-Dosprey.bench.outputDir=' + str(root),
       '-Dosprey.bench.pdbPath=' + str(pdb), '-Dosprey.bench.mutable=' + row['mutable'],
       '-Dosprey.bench.flexible=' + row['flexible'],
       '-cp', cp, 'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
(root / 'design.json').write_text(json.dumps(row, indent=2) + '\n')
shutil.copy2(prep / 'package/protocol.json', root / 'protocol.original.json')
protocol = json.loads((root / 'protocol.original.json').read_text())
protocol.update(memory_gib=128, heap_gib=96,
                resource_override='Lower memory for previously pending tasks; inputs and algorithm unchanged',
                replaces_job=f'12632191_{task}', original_prep_root=str(prep))
(root / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
(root / 'command.sh').write_text(shlex.join(cmd) + '\n')
for name, args in [('slurm_job.txt', ['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID']]),
                   ('lscpu.txt', ['lscpu']),
                   ('jvm_preflight.txt', [java, '-Xmx96g', '-Xms4g', '-XX:ActiveProcessorCount=16', '-version'])]:
    with (root / name).open('w') as f:
        subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)
manifest = dict(design_id=row['design_id'], system=row['system'], task=task,
                slurm_job_id=os.environ['SLURM_JOB_ID'], node=os.environ.get('SLURMD_NODENAME'),
                prep_root=str(prep), start_epoch=time.time(), status='RUNNING', precision='FP64',
                ccd='CPU', gpus=0, cpus=16, memory_gib=128, heap_gib=96,
                replaces_job=f'12632191_{task}',
                partition=os.environ.get('SLURM_JOB_PARTITION'), stability_filter=False, fresh_emats=True)
(root / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('MARKSTAR_START', row['design_id'], 'flex_delta', row['flex_delta'], 'root', root, flush=True)
with (root / 'run.log').open('w') as log:
    result = subprocess.run(['/usr/bin/time', '-o', str(root / 'wall.time'),
                             '-f', 'elapsed=%e maxRssKiB=%M', *cmd],
                            cwd=root, stdout=log, stderr=subprocess.STDOUT)
manifest.update(end_epoch=time.time(), exit_code=result.returncode,
                status='COMPLETED' if result.returncode == 0 else 'FAILED')
path = root / (row['design_id'] + '_markstar.csv')
if path.exists():
    with path.open() as f:
        rows = list(csv.DictReader(f))
    manifest.update(result_rows=len(rows),
                    all_estimated=sum(all(r[s + '_status'] == 'Estimated'
                                          for s in ('prot', 'lig', 'comp')) for r in rows))
    if len(rows) != row['expected_sequences']:
        manifest['status'] = 'INCOMPLETE'
else:
    manifest['status'] = 'MISSING_RESULT' if result.returncode == 0 else 'FAILED'
(root / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('MARKSTAR_END', json.dumps(manifest), flush=True)
sys.exit(result.returncode or (0 if manifest['status'] == 'COMPLETED' else 5))
