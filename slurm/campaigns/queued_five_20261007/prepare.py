"""Freeze and preflight the five approved MARK* additions inside Slurm."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm under grisman.'
HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE/'plan.json').read_text())
assert not (HERE/'registry.json').exists(), 'Inspect the existing preparation/submission before repeating.'
OUT = Path('/usr/xtmp/lz280/markstar_queued_five_20261007')/('prep_'+os.environ['SLURM_JOB_ID'])
BASE = Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
SUPPLEMENT = Path(json.loads((HERE.parent/'plan_supplements_20261007/registry.json').read_text())['package'])

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

definitions = {}
for source in (BASE, SUPPLEMENT):
    for row in read(source/'designs.tsv'):
        definitions[row['design_id']] = (row, source)
historical = {r['design']: r for r in read(Path('/usr/xtmp/lz280/markstar_extension_inventory_12826374/all161.tsv'))}
inventory = HERE.parent/'current_inventory_20261006'
table = {r['design']: r for r in read(inventory/'current_plan.tsv')}
prefs = json.loads((inventory/'reporting_preferences.json').read_text())
excluded = set(json.loads((inventory/'excluded_designs.json').read_text())['designs'])
requested = {r['design'] for r in PLAN['requests']}
assert len(requested)==5 and not requested & (excluded | set(prefs['conditional_designs']))
for name in requested:
    assert table[name]['state']=='待提交' and table[name]['recommendation']=='保留', table[name]
    assert not table[name]['job']

sources = [source/row['pdb_relative'] for name in requested for row, source in [definitions[name]]]
print('STORAGE_ESTIMATE', json.dumps(dict(pdb_files=5, pdb_bytes=sum(p.stat().st_size for p in sources),
      preparation_files_under=100, preparation_bytes_under=50*2**20,
      runtime='Five independent MARK* runs; existing classpaths reused; all matrices, logs and outputs under scratch.')), flush=True)
OUT.mkdir(parents=True, exist_ok=False)
package = OUT/'package'
(package/'structures').mkdir(parents=True)
(package/'preflight').mkdir()
shutil.copy2(__file__, OUT/'prepare.py')
shutil.copy2(HERE/'plan.json', OUT/'plan.json')
packs = {}
for group in (Path('/usr/xtmp/lz280/packstar_gpu101_20261005/run_12814879'),
              Path('/usr/xtmp/lz280/markstar_plan_supplements_20261007/pack_12832109')):
    for row in read(group/'status.tsv'):
        if row['design'] in requested and row['status'] in ('COMPLETED', 'INCOMPLETE_ESTIMATES'):
            packs[row['design']] = row
assert set(packs)==requested
configs, rows, verification = {}, [], []
for index, request in enumerate(PLAN['requests']):
    name = request['design']
    original, source_package = definitions[name]
    row = dict(original)
    source_pdb = source_package/row['pdb_relative']
    assert digest(source_pdb)==row['pdb_sha256']
    pdb = package/row['pdb_relative']
    shutil.copy2(source_pdb, pdb)
    assert digest(pdb)==row['pdb_sha256']
    anchor = name if name in historical else row['anchor_design']
    source = Path(historical[anchor]['mark_path'])
    command = shlex.split((source/'command.sh').read_text())
    updates = {'-XX:ActiveProcessorCount=':'64', '-Dosprey.bench.numCPUs=':'64',
               '-Dosprey.bench.designId=':name, '-Dosprey.bench.outputDir=':'{RUN}',
               '-Djava.io.tmpdir=':'{RUN}/tmp', '-Dosprey.bench.pdbPath=':str(pdb),
               '-Dosprey.bench.mutable=':row['mutable'], '-Dosprey.bench.flexible=':row['flexible']}
    rewritten, changed = [], set()
    for token in command:
        key = next((k for k in updates if token.startswith(k)), None)
        if key:
            token = key+updates[key]
            changed.add(key)
        rewritten.append(token.replace(str(source), '{RUN}'))
    assert changed==set(updates) and '-Dosprey.bench.method=markstar' in rewritten
    assert all(Path(p).exists() for p in rewritten[rewritten.index('-cp')+1].split(':'))
    heap = int(re.fullmatch(r'-Xmx(\d+)g', next(t for t in rewritten if t.startswith('-Xmx')))[1])
    assert 0 < heap < PLAN['memory_gib']
    pack_root = Path(packs[name]['output'])
    pack_manifest = json.loads((pack_root/'manifest.json').read_text())
    for key in ('mutable', 'flexible', 'pdb_sha256', 'expected_sequences'):
        assert str(pack_manifest['design'][key])==str(row[key]), (name, key)
    pack_results = list(csv.DictReader((pack_root/(name+'_packstar.csv')).open()))
    assert len(pack_results)==int(row['expected_sequences'])
    pack_seconds = float(re.search(r'elapsed=([\d.]+)', (pack_root/'wall.time').read_text())[1])
    assert pack_seconds < 1800, (name, pack_seconds)
    preflight = package/'preflight'/name
    preflight.mkdir()
    (preflight/'tmp').mkdir()
    precommand = []
    for token in rewritten:
        token = token.replace('{RUN}', str(preflight))
        if token.startswith('-Xmx'):
            token = '-Xmx8g'
        elif token.startswith('-XX:ActiveProcessorCount='):
            token = '-XX:ActiveProcessorCount=4'
        elif token.startswith('-Dosprey.bench.numCPUs='):
            token = '-Dosprey.bench.numCPUs=4'
        elif token.startswith('-Dosprey.bench.method='):
            token = '-Dosprey.bench.method=sequence_dump'
        precommand.append(token)
    (preflight/'command.sh').write_text(shlex.join(precommand)+'\n')
    with (preflight/'run.log').open('w') as log:
        subprocess.run(precommand, cwd=preflight, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    log = (preflight/'run.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log
    positions = re.findall(r'\[FRONTIER_POSITION\] state=Complex residue=(\S+)', log)
    assert len(positions)==int(row['total_positions'])
    assert set(positions)==set((row['mutable']+';'+row['flexible']).split(';'))
    sequences = read(preflight/(name+'_sequences.tsv'))
    assert len(sequences)==int(row['expected_sequences'])
    assert {r['sequence'] for r in sequences}=={r['sequence'] for r in pack_results}
    row.update(destination=request['group'], group_index=index, requested_cpus=64, requested_memory_gib=192,
               java_heap_gib=heap, limit_hours=336, kind='USER_APPROVED_PLAN_ADDITION',
               reference_design=anchor, reference_mark_path=str(source), source_package=str(source_package))
    rows.append(row)
    configs[name] = dict(row=row, command_template=rewritten, input_pdb=str(pdb),
        model_contains='AMD EPYC 9554' if request['group']=='compsci' else 'Xeon(R) Gold 5320',
        source_command_sha256=digest(source/'command.sh'))
    verification.append(dict(design=name, preflight='PASS', sequence_set_matches_PACK=True,
        expected_sequences=len(sequences), positions=len(positions), pack_seconds=pack_seconds,
        pack_estimated=int(packs[name]['estimated_rows']), pack_output=str(pack_root), heap_gib=heap))
    print('PREFLIGHT_OK', json.dumps(verification[-1]), flush=True)
write(package/'designs.tsv', rows)
save(package/'run_configs.json', configs)
for filename in ('markstar_runner.py', 'markstar_array.slurm'):
    if filename.endswith('.py'):
        compile((HERE/filename).read_text(), filename, 'exec')
    else:
        subprocess.run(['bash', '-n', str(HERE/filename)], check=True)
    shutil.copy2(HERE/filename, OUT/filename)
(OUT/'READY').write_text('Five MARK* input preflights passed and sequence sets match measured PACK* results.\n')
write(HERE/'designs.tsv', rows)
record = dict(status='PREPARED', management_job=os.environ['SLURM_JOB_ID'], prep=str(OUT),
              package=str(package), verification=verification, requests=PLAN['requests'], jobs=[])
save(HERE/'registry.json', record)
save(OUT/'preparation.json', record)
print('PREPARED', OUT, flush=True)
