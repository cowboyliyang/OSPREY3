"""Audit active input pairs, append-only manifests and scheduler dependencies."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here = Path(__file__).parent
repo = here.parents[2]
previous = here.parent/'deadline_focus_20261006'
registry = json.loads((here/'registry.json').read_text())
prep = Path(registry['prep_root'])

def read(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))

def original(path):
    text = subprocess.check_output(['git','show','e3caa95f:'+str(path.relative_to(repo))], cwd=repo, text=True)
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))

rows = read(prep/'package/designs.tsv')
assert rows == read(here/'pack_designs.tsv') == read(here/'new_designs.tsv')
assert read(here/'requested_designs.tsv') == read(prep/'package/provenance/requested_designs.tsv')
assert len(rows)==9 and all(r['pack_platform']=='gpu' for r in rows)
assert hashlib.sha256((prep/'package/designs.tsv').read_bytes()).hexdigest()==registry['designs_sha256']
for root in [prep, Path('/usr/xtmp/lz280/markstar_near14_20261006/prep_12827301'),
             Path('/usr/xtmp/lz280/markstar_near14_20261006/prep_12827308')]:
    subprocess.run(['sha256sum','-c','SHA256SUMS'], cwd=root/'package', stdout=subprocess.DEVNULL, check=True)
for name in ['run_group.py','run_group.slurm']:
    assert (here/name).read_bytes() == (prep/name).read_bytes()
assert not (here/'run_cpu_group.slurm').exists()
by_name = {r['design_id']:r for r in rows}
assert '2p4a_flex_p18' not in by_name
r = by_name['2p4a_flex_p12']
assert r['expected_sequences']=='39' and r['historical_scan_branchwidth']=='5'
assert r['added_since_anchor'].split(';')==['C356','D469','C357']
assert r['pack_reference_platform']=='cpu' and by_name['3bu8_flex_p12']['pack_reference_platform']=='cpu'
protocol = json.loads((prep/'package/protocol.json').read_text())
assert protocol['platforms']=={'gpu':dict(cpus=128,gpus=4,model='RTX PRO 6000',global_limit_seconds=21600)}
properties = (prep/'build/source/slurm/h200/production.properties').read_text()
assert 'branchdp.dp.gpu=true' in properties and 'branchdp.dp.gpu.failIfNoGpuPath=true' in properties
pack_detail = subprocess.check_output(['scontrol','show','job',registry['packstar_job'],'-o'], text=True)
assert 'CPUs/Task=128 ' in pack_detail and 'gres/gpu:rtx_pro_6000=4' in pack_detail
assert 'TimeLimit=06:00:00 ' in pack_detail and 'afterany:12827169' in pack_detail

markprep = Path(registry['markstar_prep'])
configs = json.loads((markprep/'package/run_configs.json').read_text())
assert set(configs)==set(by_name)
for name, config in configs.items():
    command = config['command_template']
    row = by_name[name]
    assert '-Dosprey.bench.numCPUs=64' in command and '-XX:ActiveProcessorCount=64' in command
    for key in ['mutable','flexible']:
        assert '-Dosprey.bench.'+key+'='+row[key] in command
    assert hashlib.sha256(Path(config['input_pdb']).read_bytes()).hexdigest()==row['pdb_sha256']
    if row['destination']=='compsci':
        assert config['model_contains']=='AMD EPYC 9554'
assert 'FOCUS_PACK_ROOT' not in (markprep/'markstar_runner.py').read_text()
for job in registry['local_markstar_jobs']:
    detail = subprocess.check_output(['scontrol','show','job',job['job'],'-o'], text=True)
    assert 'CPUs/Task=64 ' in detail and 'afterany:'+job['after_job'] in detail
    assert 'TimeLimit=14-00:00:00 ' in detail and 'Account=grisman ' in detail

old_edges = original(previous/'local_lane_dependencies.tsv')
edges = read(previous/'local_lane_dependencies.tsv')
assert edges[:len(old_edges)]==old_edges and len(edges)==len(old_edges)+4
assert len({r['parent'] for r in edges})==len(edges), 'No lane may branch into concurrent successors'
for r in old_edges:
    detail = subprocess.check_output(['scontrol','show','job',r['job'],'-o'], text=True)
    if 'JobState=PENDING' in detail:
        assert 'afterany:'+r['parent'] in detail, (r,detail)
old_priority = original(previous/'dcc_launch_priority.tsv')
priority = read(previous/'dcc_launch_priority.tsv')
assert priority[:14]==old_priority and len(priority)==19
assert [r['design_id'] for r in priority[-5:]]==[r['design_id'] for r in rows if r['destination']=='DCC']
dcc_new = read(previous/'dcc_new_designs.tsv')
old_dcc = original(previous/'dcc_new_designs.tsv')
assert len(dcc_new)==9 and len({r['design_id'] for r in dcc_new})==9
for old, current in zip(old_dcc,dcc_new):
    assert all(current[key]==value for key,value in old.items())
old_actions = original(previous/'design_actions.tsv')
actions = read(previous/'design_actions.tsv')
assert len(actions)==len(old_actions)+9 and len({r['design'] for r in actions})==len(actions)
for old, current in zip(old_actions,actions):
    if old['design']!='2xxm_flex_p1':
        assert old==current
assert not any(r['design']=='2p4a_flex_p18' for r in actions)
external = json.loads((previous/'external_status_snapshot.json').read_text())
assert external['completed_seconds']['2xxm_flex_p1']==6338 and external['reported_completed_count']==30
assert registry['existing_jobs_modified']==[] and registry['dcc_status']=='NOT_APPLIED_REMOTELY'
for path in here.glob('*.py'):
    compile(path.read_text(), str(path), 'exec')
for path in here.glob('*.slurm'):
    subprocess.run(['bash','-n',str(path)], check=True)
subprocess.run(['git','diff','--check'], cwd=repo, check=True)

history = json.loads(Path('/usr/xtmp/lz280/markstar_extension_inventory_12826374/inventory.json').read_text())['systems']['2p4a']
summary = dict(verification_job=os.environ['SLURM_JOB_ID'], packstar_job=registry['packstar_job'],
    gpu_designs=9, local_markstar_jobs=4, dcc_additions_handoff_only=5,
    original_queue_and_definitions_preserved=True, old_frozen_packages_verified=True,
    pack_and_mark_inputs_match=True, gpu_request_and_dependency_verified=True,
    conservative_design='2p4a_flex_p12', additional_residues=3, expected_sequences=39,
    two_p4a_history=[{k:r[k] for k in ['design','kind','mark_hours','pack_minutes','ratio','historical_cpu']} for r in history])
(here/'verification.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2), flush=True)
