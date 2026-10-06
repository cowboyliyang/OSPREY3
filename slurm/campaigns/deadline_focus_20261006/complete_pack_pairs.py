"""Pair restored CPU continuation tiers with PACK* in the pending GPU group.

Create a fresh package; preserve all original frozen inputs and CPU jobs.
If either old GPU job has started, stop before scheduler changes.
"""
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm'
here = Path(__file__).parent
registry = json.loads((here/'registry.json').read_text())
assert 'pack_pair_completion' not in registry, 'Already applied; inspect the registry'
original = Path(registry['prep_root'])
extension = Path('/usr/xtmp/lz280/frontier_extension7_20261006/prep_12826503')
root = original.parent/('pack_pairs_'+os.environ['SLURM_JOB_ID'])
root.mkdir(exist_ok=False)
package = root/'package'
additions = ['2rl0_flex_p4', '4u3s_flex_p2']
old_jobs = [registry['packstar_job'], registry['packstar_screening_job']]


def read(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))


def write(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def detail(job):
    return subprocess.check_output(['scontrol', 'show', 'job', job, '-o'], text=True)


def queue():
    return subprocess.check_output(['squeue', '-r', '-h', '-u', 'lz280', '-o', '%i|%T|%j|%E|%C'], text=True)


for job in old_jobs:
    info = detail(job)
    assert 'JobState=PENDING ' in info and 'RunTime=00:00:00 ' in info and 'Restarts=0 ' in info, info
    (root/(job+'_before.txt')).write_text(info)
(root/'queue_before.txt').write_text(queue())

# Metadata-only storage estimate before copying small prepared inputs.
files = [p for p in (original/'package').rglob('*') if p.is_file()]
estimate = dict(package_files=len(files), package_bytes=sum(p.stat().st_size for p in files),
                extra_preflight_bytes_under=10*2**20, extra_preflight_files_under=50,
                production='Two additional workloads under xtmp; ordinary logs/caches expected in GiB. Reuse compiled classes without copying them.')
(root/'storage_estimate.json').write_text(json.dumps(estimate, indent=2)+'\n')
print('STORAGE_ESTIMATE', json.dumps(estimate), flush=True)
for source in [original, extension]:
    assert (source/'READY').exists()
    with (root/(source.name+'_checksum.log')).open('w') as f:
        subprocess.run(['sha256sum', '-c', 'SHA256SUMS'], cwd=source/'package', stdout=f, stderr=subprocess.STDOUT, check=True)
for relative in ['build/source.sha256', 'build/source/slurm/h200/production.properties', 'portable_run.py']:
    assert digest(original/relative) == digest(extension/relative), relative
shutil.copytree(original/'package', package)
shutil.copytree(original/'build', root/'build')
shutil.copy2(original/'portable_run.py', root/'portable_run.py')
for name in ['run_group.py', 'run_group.slurm', 'screen_pack.py', 'screen_pack.slurm']:
    source = here/name
    if name.endswith('.py'):
        compile(source.read_text(), str(source), 'exec')
    else:
        subprocess.run(['bash', '-n', str(source)], check=True)
    shutil.copy2(source, root/name)
shutil.copy2(__file__, root/Path(__file__).name)

rows = read(package/'designs.tsv')
assert len(rows) == 10 and not set(additions) & {r['design_id'] for r in rows}
original_rows = [dict(r) for r in rows]
available = {r['design_id']: r for r in read(extension/'package/designs.tsv')}
cpu = json.loads((Path(registry['markstar_prep'])/'package/run_configs.json').read_text())
verification = json.loads((original/'verification.json').read_text())['preflight']
previous_checks = {r['design_id']: r for r in json.loads((extension/'verification.json').read_text())['preflight']}
for name in additions:
    row = dict(available[name], task_id=str(len(rows)), destination=cpu[name]['row']['destination'], kind='retained_extension')
    for key in ['mutable', 'flexible', 'pdb_sha256', 'expected_sequences', 'total_positions']:
        assert str(row[key]) == str(cpu[name]['row'][key]), (name, key)
    for key in ['mutable', 'flexible']:
        assert '-Dosprey.bench.'+key+'='+row[key] in cpu[name]['command_template']
    assert digest(package/row['pdb_relative']) == digest(Path(cpu[name]['input_pdb'])) == row['pdb_sha256']
    assert previous_checks[name]['mark_preflight'] == previous_checks[name]['pack_preflight'] == 'PASS'
    shutil.copytree(extension/'package/preflight'/name, package/'preflight'/name)
    rows.append(row)
assert rows[:10] == original_rows and len(rows) == 12
write(package/'designs.tsv', rows)
(package/'designs.json').write_text(json.dumps(rows, indent=2)+'\n')
requests = read(package/'provenance/requested_designs.tsv')
requests += [{key: next(r for r in rows if r['design_id']==name)[key] for key in requests[0]} for name in additions]
write(package/'provenance/requested_designs.tsv', requests)

# Re-run both added PACK* input preflights against the exact current build.
env = dict(os.environ, JAVA_HOME='/home/users/lz280/java/jdk-17.0.2+8', BUILD_ROOT=str(root/'build'),
           INPUT_ROOT=str(package), RESULT_ROOT=str(root/'pack_preflight'), PACKSTAR_EXPECT_GPU='RTX PRO 6000')
for row in rows[10:]:
    name = row['design_id']
    with (root/(name+'_preflight.log')).open('w') as f:
        subprocess.run(['python3', str(root/'portable_run.py'), '--frontier-manifest', str(package/'designs.tsv'),
                        '--index', row['task_id'], '--mode', 'preflight', '--arm', 'pair-only', '--rb', '1',
                        '--seed', '42', '--gpus', '4', '--gpu-gib', '85', '--heap-gib', '8', '--host-gib', '4'],
                       env=env, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=240)
    manifests = list((root/'pack_preflight').glob(name+'_preflight_*/manifest.json'))
    assert len(manifests)==1
    manifest = json.loads(manifests[0].read_text())
    assert manifest['status']=='COMPLETED' and manifest['result_rows']==int(row['expected_sequences'])
    result = manifests[0].parent
    log = (result/'run.log').read_text()
    assert int(re.search(r'Complex positions: (\d+)', log)[1])==int(row['total_positions'])
    pack_sequences = read(result/(name+'_sequences.tsv'))
    mark_sequences = read(package/'preflight'/name/(name+'_sequences.tsv'))
    assert [r['sequence'] for r in pack_sequences] == [r['sequence'] for r in mark_sequences]
    verification.append(dict(previous_checks[name], paired_cpu_command_matches=True, pack_preflight='PASS'))
    print('PAIR_PREFLIGHT_OK', name, len(pack_sequences), flush=True)

protocol = json.loads((package/'protocol.json').read_text())
protocol.update(designs=12, systems=len({r['system'] for r in rows}),
                scheduling='One exclusive node allocation; twelve sequential workloads; six hours total',
                markstar_status='CPU jobs are independent of PACK* measurements')
assert protocol['global_limit_seconds']==21600 and protocol['gpus']==4
assert protocol['global_limit_seconds']-protocol['shutdown_reserve_seconds'] >= len(rows)*protocol['future_case_reserve_seconds']
(package/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
(root/'verification.json').write_text(json.dumps(dict(preflight=verification, protocol=protocol,
    designs_sha256=digest(package/'designs.tsv')), indent=2)+'\n')
checks = [digest(p)+'  '+str(p.relative_to(package)) for p in sorted(package.rglob('*')) if p.is_file() and p.name!='SHA256SUMS']
(package/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(root/'READY').write_text('Twelve PACK* designs; two restored continuation pairs verified against active CPU commands\n')

# Swap only never-started GPU jobs; validate replacements while held.
held = []
submitted = []
retired = False
try:
    for job in old_jobs:
        assert 'JobState=PENDING ' in detail(job)
        subprocess.run(['scontrol', 'hold', job], check=True)
        held.append(job)
        info = detail(job)
        assert all(word in info for word in ['JobState=PENDING ', 'RunTime=00:00:00 ', 'Restarts=0 ', 'Reason=JobHeldUser ']), info
    os.environ['FOCUS_PACK_INDICES'] = ','.join(str(i) for i in range(12))
    def submit(args):
        job = subprocess.check_output(['sbatch', '--parsable', '--hold', '--account=grisman', *args], text=True).strip().split(';')[0]
        assert job.isdigit()
        submitted.append(job)
        (root/'submitted_jobs.json').write_text(json.dumps(submitted)+'\n')
        return job
    pack_job = submit(['--export=ALL,EXTENSION_PREP='+str(root), str(root/'run_group.slurm')])
    pack_output = root.parent/('pack_'+pack_job)
    screen_job = submit(['--dependency=afterany:'+pack_job,
        '--export=ALL,EXTENSION_PREP='+str(root)+',FOCUS_PACK_ROOT='+str(pack_output)+',FOCUS_PACK_JOB='+pack_job,
        str(root/'screen_pack.slurm')])
    info = detail(pack_job)
    assert all(word in info for word in ['JobState=PENDING ', 'Account=grisman ', 'CPUs/Task=128 ',
               'TimeLimit=06:00:00 ', 'Partition=compsci-gpu ', 'gres/gpu:rtx_pro_6000=4']), info
    assert 'afterany:'+pack_job in detail(screen_job)
    for job in old_jobs:
        assert 'JobState=PENDING ' in detail(job)
    for job in reversed(old_jobs):
        subprocess.run(['scancel', '--state=PENDING', job], check=True)
    retired = True
    for job in old_jobs:
        assert 'JobState=CANCELLED ' in detail(job), job
except BaseException:
    if not retired:
        for job in submitted:
            subprocess.run(['scancel', '--state=PENDING', job], check=False)
        for job in held:
            subprocess.run(['scontrol', 'release', job], check=False)
    raise

record = dict(applied_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), audit=str(root),
              added_designs=additions, original_prep=str(original), retired_pending_jobs=old_jobs,
              pack_job=pack_job, screen_job=screen_job, designs=12, gpus=4, cpus=128, limit_hours=6,
              first_ten_definitions_and_order_preserved=True, paired_cpu_inputs_and_sequences_verified=True,
              cpu_jobs_unchanged=True, status='SUBMITTED_HELD')
registry.update(prep_job=os.environ['SLURM_JOB_ID'], prep_root=str(root), package=str(package), design_count=12,
                designs_sha256=digest(package/'designs.tsv'), verification=verification,
                packstar_job=pack_job, packstar_replaced_pending_job=old_jobs[0],
                packstar_output=str(pack_output), packstar_screening_job=screen_job,
                packstar_selected_indices=os.environ['FOCUS_PACK_INDICES'], pack_pair_completion=record)
write(here/'pack_designs.tsv', rows)
write(here/'requested_designs.tsv', requests)
(here/'registry.json').write_text(json.dumps(registry, indent=2)+'\n')
(here/'pack_pair_completion.json').write_text(json.dumps(record, indent=2)+'\n')
for job in [screen_job, pack_job]:
    subprocess.run(['scontrol', 'release', job], check=True)
record['status'] = 'SUBMITTED'
(here/'registry.json').write_text(json.dumps(registry, indent=2)+'\n')
(here/'pack_pair_completion.json').write_text(json.dumps(record, indent=2)+'\n')
(root/'queue_after.txt').write_text(queue())
print(json.dumps(record, indent=2), flush=True)
