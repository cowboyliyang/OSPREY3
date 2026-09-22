"""Publish the active H200 manifest and a checked portable input package on Slurm."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

assert os.environ.get('SLURM_JOB_ID'), 'Submit through Slurm'
prep, root = map(Path, sys.argv[1:3])
array = sys.argv[3]
assert array.isdigit() and (prep / 'READY').read_text().strip() == '25 designs validated'
repo = Path('/home/users/lz280/OSPREY3-fresh-packstar')
config = repo / 'slurm/h200'
package = root / 'package'
package.mkdir()
for name in ['structures', 'preflight', 'provenance']:
    (package / name).mkdir()


def read(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))


def write(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def verify(source):
    for line in (source / 'SHA256SUMS').read_text().splitlines():
        checksum, relative = line.split('  ', 1)
        assert hashlib.sha256((source / relative).read_bytes()).hexdigest() == checksum, (source, relative)


sources = {
    'frontier.tsv': Path('/usr/xtmp/lz280/packstar_flex_frontier12_20260918/handoff_12632294/package'),
    'frontier_add4.tsv': Path('/usr/xtmp/lz280/packstar_flex_frontier16_20260918/handoff_12633791/package'),
    'frontier_small4.tsv': Path('/usr/xtmp/lz280/packstar_flex_frontier20_20260919/handoff_12643037/package'),
    'frontier_remaining18_plus4.tsv': Path('/usr/xtmp/lz280/packstar_flex_frontier38_20260920/prep_12651495/package'),
    'frontier_backfill25.tsv': prep / 'package',
}
for source in sources.values():
    verify(source)
fields = list(read(config / 'frontier.tsv')[0])
added = []
for design in json.loads((prep / 'package/designs.json').read_text()):
    row = {k: design.get(k, '') for k in fields}
    row.update(markstar_job=f'{array}_{design["task_id"]}', source_package=str(prep / 'package'),
               source_task_id=design['task_id'])
    added.append(row)
write(root / 'frontier_backfill25.tsv', added)
excluded = {'3k3q_flex_p1', '3k3q_flex_p2'}
rows = []
for name, source in sources.items():
    path = root / name if name == 'frontier_backfill25.tsv' else config / name
    batch = read(path)
    source_rows = {r['design_id']: r for r in read(source / 'designs.tsv')}
    for r in batch:
        frozen = source_rows[r['design_id']]
        for key in ['mutable', 'flexible', 'pdb_sha256', 'expected_sequences', 'total_positions']:
            assert r[key] == frozen[key], (r['design_id'], key)
        if r['design_id'] in excluded:
            continue
        rows.append(dict(r, task_id=len(rows), h200_source_manifest=name, h200_source_index=r['task_id']))
assert len(rows) == 162 and len({r['design_id'] for r in rows}) == 162
assert len({r['system'] for r in rows}) == 38
assert [r['design_id'] for r in rows if r['system'] == '3k3q'] == ['3k3q_flex_p0']
assert all(r['h200_source_manifest'] == 'frontier_backfill25.tsv' for r in rows[137:])
pdb_paths = {r['pdb_relative']: sources[r['h200_source_manifest']] / r['pdb_relative'] for r in rows}
print('STORAGE_ESTIMATE', json.dumps(dict(pdb_files=len(pdb_paths), pdb_bytes=sum(p.stat().st_size for p in pdb_paths.values()),
      designs=162, preflight_files=324, note='Small portable input package; no builds, EMAT or run caches.')), flush=True)
seqsets = {}
for row in rows:
    source = sources[row['h200_source_manifest']]
    pdb = source / row['pdb_relative']
    assert hashlib.sha256(pdb.read_bytes()).hexdigest() == row['pdb_sha256']
    dest = package / row['pdb_relative']
    if not dest.exists():
        shutil.copy2(pdb, dest)
    assert hashlib.sha256(dest.read_bytes()).hexdigest() == row['pdb_sha256']
    preflight = source / 'preflight' / row['design_id']
    seqname = row['design_id'] + '_sequences.tsv'
    seqs = read(preflight / seqname)
    assert len(seqs) == int(row['expected_sequences'])
    values = [r['sequence'] for r in seqs]
    if row['system'] in seqsets:
        assert values == seqsets[row['system']], row['design_id']
    seqsets[row['system']] = values
    log = (preflight / 'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log
    positions = [x.split(' residue=')[1].split()[0] for x in log.splitlines() if x.startswith('[FRONTIER_POSITION] state=Complex ')]
    assert len(positions) == int(row['total_positions']) and set(positions) == set((row['mutable'] + ';' + row['flexible']).split(';'))
    out = package / 'preflight' / row['design_id']
    out.mkdir()
    for name in [seqname, 'preflight.log']:
        shutil.copy2(preflight / name, out / name)
write(package / 'designs.tsv', rows)
write(root / 'frontier_active.tsv', rows)
protocol = dict(designs=162, systems=38, source_packages={k: str(v) for k, v in sources.items()},
    excluded_designs=sorted(excluded), retained_3k3q='3k3q_flex_p0', added_markstar_array=array,
    h200_index='Independent contiguous 0..161; align old results by design_id.',
    h200_resources='Use the already validated H200 PACK* build and its uniform resource/protocol settings.',
    markstar_backfill_resources=dict(cpus=16, memory_gib=192, heap_gib=160, partition='compsci', account='grisman', days=14))
(package / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
(package / 'README.md').write_text('# Current 38-system H200 inputs\n\n162 designs, index 0..161. 3K3Q +1/+2 removed; 0 retained.\nIndices 137..161 are the new 25 extensions. Use run_active.slurm with the matching frontier_active.tsv.\nMatch existing results by design_id before selecting missing tasks.\n')
# Validate the actual new launcher using the committed portable runner and an
# existing compiled CPU protocol. Keep unrelated working-tree runner edits out.
validation_repo = root / 'validation_repo'
vconfig = validation_repo / 'slurm/h200'
vconfig.mkdir(parents=True)
(vconfig / 'run.py').write_bytes(subprocess.check_output(['git', '-C', str(repo), 'show', 'HEAD:slurm/h200/run.py']))
shutil.copy2(config / 'run_active.slurm', vconfig / 'run_active.slurm')
shutil.copy2(root / 'frontier_active.tsv', vconfig / 'frontier_active.tsv')
subprocess.run(['bash', '-n', str(vconfig / 'run_active.slurm')], check=True)
build = root / 'validation_build'
(build / 'source/slurm/h200').mkdir(parents=True)
historical = Path('/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626762')
for name in ['git_head.txt', 'source.sha256']:
    shutil.copy2(historical / name, build / name)
shutil.copy2(prep / 'test_classpath.txt', build / 'test_classpath.txt')
shutil.copy2(config / 'production.properties', build / 'source/slurm/h200/production.properties')
(build / 'READY').write_text('CPU validation only\n')
selected = [r for r in rows if int(r['task_id']) >= 137]
selected += [rows[0], next(r for r in rows if r['system'] == '3cal'), next(r for r in rows if r['system'] == '4wem')]
env = dict(os.environ, REPO=str(validation_repo), BUILD_ROOT=str(build), INPUT_ROOT=str(package),
           RESULT_ROOT=str(root / 'runner_validation'), JAVA_HOME='/home/users/lz280/java/jdk-17.0.2+8')
for row in selected:
    with (root / f'runner_{row["task_id"]}.log').open('w') as log:
        subprocess.run(['bash', str(vconfig / 'run_active.slurm'), '--mode', 'preflight', '--heap-gib', '8', '--host-gib', '4'],
            env=dict(env, SLURM_ARRAY_TASK_ID=str(row['task_id'])), stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    print('H200_LAUNCHER_PREFLIGHT_OK', row['task_id'], row['design_id'], flush=True)
checks = [hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + str(p.relative_to(package)) for p in sorted(package.rglob('*')) if p.is_file()]
(package / 'SHA256SUMS').write_text('\n'.join(checks) + '\n')
audit = dict(protocol, prep_root=str(prep), handoff_root=str(root), validation_job=os.environ['SLURM_JOB_ID'],
    validation='All 162 frozen position lists, PDB checksums and sequence lists verified; 28 actual CPU preflights through run_active.slurm, including all 25 additions, 3CAL and 4WEM. No H200 GPU execution.',
    package_files=sum(p.is_file() for p in package.rglob('*')),
    package_bytes=sum(p.stat().st_size for p in package.rglob('*') if p.is_file()),
    sequences=sum(int(r['expected_sequences']) for r in rows),
    active_manifest_sha256=hashlib.sha256((root / 'frontier_active.tsv').read_bytes()).hexdigest())
(root / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
# Only small versioned manifests and audit summaries are published in home.
for name in ['frontier_backfill25.tsv', 'frontier_active.tsv']:
    target = config / name
    assert not target.exists(), target
    shutil.copy2(root / name, target)
target = config / 'frontier_active.audit.json'
assert not target.exists()
shutil.copy2(root / 'audit.json', target)
(root / 'READY').write_text('162 active designs ready for H200 transfer\n')
print('READY', json.dumps(audit), flush=True)
