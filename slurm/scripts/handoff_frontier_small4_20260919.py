"""Create and validate the independent H200 add-on package, through Slurm."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

assert os.environ.get('SLURM_JOB_ID')
prep, root = map(Path, sys.argv[1:3])
small_job, large_job = sys.argv[3:5]
repo = Path('/home/users/lz280/OSPREY3-fresh-packstar')
assert (prep / 'READY').read_text().strip() == '16 designs validated'
source = prep / 'package'
for line in (source / 'SHA256SUMS').read_text().splitlines():
    checksum, relative = line.split('  ', 1)
    assert hashlib.sha256((source / relative).read_bytes()).hexdigest() == checksum
root.mkdir(parents=True, exist_ok=False)
package = root / 'package'
package.mkdir()
for name in ('structures', 'preflight', 'config'):
    shutil.copytree(source / name, package / name)
for name in ('protocol.json', 'README.md'):
    shutil.copy2(source / name, package / name)
original_manifest = repo / 'slurm/h200/frontier.tsv'
fields = next(csv.reader(original_manifest.open(), delimiter='\t'))
designs = json.loads((source / 'designs.json').read_text())
rows = []
for design in designs:
    row = {key: design.get(key, '') for key in fields}
    task = design['task_id']
    row.update(markstar_job=f'{large_job if design["system"] == "3bua" else small_job}_{task}',
               markstar_cpus=16, source_package=str(source), source_task_id=task)
    rows.append(row)
assert len(rows) == 16 and len({r['design_id'] for r in rows}) == 16
assert not {r['design_id'] for r in rows} & {r['design_id'] for r in csv.DictReader(original_manifest.open(), delimiter='\t')}
with (package / 'designs.tsv').open('w') as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
# Use the existing compiled Java protocol for a real CPU-only validation of
# the portable runner's independent manifest path. This is not an H200 test.
build = root / 'validation_build'
(build / 'source/slurm/h200').mkdir(parents=True)
historical_build = Path('/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626762')
for name in ('git_head.txt', 'source.sha256'):
    shutil.copy2(historical_build / name, build / name)
shutil.copy2(prep / 'test_classpath.txt', build / 'test_classpath.txt')
shutil.copy2(repo / 'slurm/h200/production.properties', build / 'source/slurm/h200/production.properties')
shutil.copy2(original_manifest, build / 'source/slurm/h200/frontier.tsv')
shutil.copy2(repo / 'slurm/h200/run.py', root / 'run.py')
for name in ('run_small4.slurm',):
    subprocess.run(['bash', '-n', str(repo / 'slurm/h200' / name)], check=True)
for name in ('prepare_flex_frontier_small4_20260919.slurm', 'run_markstar_flex_small4_20260919.slurm'):
    subprocess.run(['bash', '-n', str(repo / 'slurm/scripts' / name)], check=True)
env = dict(os.environ, BUILD_ROOT=str(build), INPUT_ROOT=str(package),
           RESULT_ROOT=str(root / 'runner_validation'),
           JAVA_HOME='/home/users/lz280/java/jdk-17.0.2+8')
for row in rows:
    cmd = [sys.executable, str(root / 'run.py'), '--frontier-manifest', str(package / 'designs.tsv'),
           '--mode', 'preflight', '--index', str(row['task_id']), '--heap-gib', '8', '--host-gib', '4']
    with (root / f'runner_{row["task_id"]}.log').open('w') as out:
        subprocess.run(cmd, env=env, stdout=out, stderr=subprocess.STDOUT, check=True, timeout=180)
    print('H200_DRIVER_PREFLIGHT_OK', row['design_id'], flush=True)
# Check the pre-existing cohort path as well, without modifying its manifest.
old_env = dict(env, INPUT_ROOT='/usr/xtmp/lz280/packstar_flex_frontier12_20260918/handoff_12632294/package',
               RESULT_ROOT=str(root / 'original_cohort_validation'))
with (root / 'original_runner.log').open('w') as out:
    subprocess.run([sys.executable, str(root / 'run.py'), '--mode', 'preflight', '--index', '0',
                    '--heap-gib', '8', '--host-gib', '4'], env=old_env,
                   stdout=out, stderr=subprocess.STDOUT, check=True, timeout=180)
checks = []
for path in sorted(package.rglob('*')):
    if path.is_file():
        checks.append(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + str(path.relative_to(package)))
(package / 'SHA256SUMS').write_text('\n'.join(checks) + '\n')
summary = dict(designs=16, systems=4, sequences=sum(int(r['expected_sequences']) for r in rows),
               pdb_bytes=sum(p.stat().st_size for p in (package / 'structures').iterdir()),
               package_files=sum(p.is_file() for p in package.rglob('*')),
               package_bytes=sum(p.stat().st_size for p in package.rglob('*') if p.is_file()),
               original_manifest_sha256=hashlib.sha256(original_manifest.read_bytes()).hexdigest(),
               validation='16 added and 1 original real sequence_dump runs via portable CPU preflight; no H200 GPU execution',
               markstar_arrays=[small_job, large_job])
(root / 'audit.json').write_text(json.dumps(summary, indent=2) + '\n')
(root / 'READY').write_text('16 additional designs ready for transfer; original cohort retained\n')
print(json.dumps(summary), flush=True)
