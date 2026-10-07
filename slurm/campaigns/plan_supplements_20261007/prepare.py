"""Freeze three missing PACK* inputs and preflight both methods through Slurm."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Submit preparation through Slurm.'
HERE = Path(__file__).resolve().parent
ROOT = Path('/usr/xtmp/lz280/markstar_plan_supplements_20261007') / ('prep_' + os.environ['SLURM_JOB_ID'])
ROOT.mkdir(parents=True, exist_ok=False)
PACKAGE = ROOT / 'package'
BASE = Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
VALIDATED = Path('/usr/xtmp/lz280/packstar_gpu101_20261005/run_12814879')
SCANS = Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config/expansion_scans')
JAVA = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

requested = read(HERE/'requested_designs.tsv')
prior = {r['design_id']: r for r in read(BASE/'designs.tsv')}
review = read(Path('/usr/xtmp/lz280/markstar_plan_review_12832102/unified_plan.tsv'))
missing = {r['design'] for r in review if not r.get('PACK_seconds')}
assert {r['design_id'] for r in requested} == missing | {'5dc0_flex_p9'}
assert not set(prior) & {r['design_id'] for r in requested}
structures = {r['system']: BASE/prior[r['anchor_design']]['pdb_relative'] for r in requested}
estimate = dict(pdb_files=len(structures), pdb_bytes=sum(p.stat().st_size for p in structures.values()),
                preparation_files_under=100, preparation_bytes_under=50*2**20,
                production='Three sequential designs; existing per-design mapped table budget 512 GiB; scratch only.')
print('STORAGE_ESTIMATE', json.dumps(estimate), flush=True)
(ROOT/'storage_estimate.json').write_text(json.dumps(estimate, indent=2)+'\n')
for name in ('structures', 'config', 'preflight', 'provenance'):
    (PACKAGE/name).mkdir(parents=True)
for system, source in structures.items():
    shutil.copy2(source, PACKAGE/'structures'/source.name)
    shutil.copy2(SCANS/(system+'_scan.tsv'), PACKAGE/'config'/(system+'_scan.tsv'))
shutil.copy2(HERE/'requested_designs.tsv', PACKAGE/'provenance/requested_designs.tsv')
shutil.copy2(__file__, ROOT/'prepare.py')
shutil.copy2(HERE.parent/'deadline_focus_20261006/run_group.py', ROOT/'run_group.py')
shutil.copy2(HERE/'run_group.slurm', ROOT/'run_group.slurm')
compile((ROOT/'run_group.py').read_text(), str(ROOT/'run_group.py'), 'exec')
subprocess.run(['bash', '-n', str(ROOT/'run_group.slurm')], check=True)
build = ROOT/'build'
(build/'source/slurm/h200').mkdir(parents=True)
for name in ('test_classpath.txt', 'git_head.txt', 'source.sha256'):
    shutil.copy2(VALIDATED/'build'/name, build/name)
shutil.copy2(VALIDATED/'build/source/slurm/h200/production.properties', build/'source/slurm/h200/production.properties')
shutil.copy2(VALIDATED/'run.py', ROOT/'portable_run.py')
assert all(Path(p).exists() for p in (build/'test_classpath.txt').read_text().strip().split(':'))

designs = []
for index, request in enumerate(requested):
    anchor = prior[request['anchor_design']]
    delta, previous_delta = int(request['flex_delta']), int(anchor['flex_delta'])
    scan = {int(r['k']): r for r in read(PACKAGE/'config'/(request['system']+'_scan.tsv'))}
    mutable = anchor['mutable'].split(';')
    anchor_flex = anchor['flexible'].split(';')
    if delta < 0:
        assert request['design_id']=='3k3q_flex_m1' and previous_delta==0
        # Same deterministic shrink rule as the existing 3bua negative tiers.
        chains = {r[0] for r in mutable}
        removed = [r for r in reversed(anchor_flex) if r[0] in chains][:1]
        assert removed == ['C439']
        flexible = [r for r in anchor_flex if r not in removed]
        all_added, new_added = [], []
    else:
        old_added = [r for r in scan[previous_delta]['added_residues'].split(';') if r]
        all_added = scan[delta]['added_residues'].split(';')
        assert all_added[:len(old_added)] == old_added
        new_added = all_added[len(old_added):]
        assert len(new_added)==delta-previous_delta
        flexible, removed = anchor_flex+new_added, []
    assert len(mutable+flexible)==len(set(mutable+flexible))==int(anchor['total_positions'])+delta-previous_delta
    pdb = PACKAGE/anchor['pdb_relative']
    assert digest(pdb)==anchor['pdb_sha256']
    residues = {line[21].strip()+line[22:26].strip()+line[26].strip()
                for line in pdb.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))}
    assert set(mutable+flexible) <= residues
    row = dict(task_id=index, design_id=request['design_id'], system=request['system'], flex_delta=delta,
               mutable=anchor['mutable'], flexible=';'.join(flexible), added=';'.join(all_added), removed=';'.join(removed),
               wt_flex_count=len(flexible), total_positions=len(mutable+flexible), expected_sequences=int(anchor['expected_sequences']),
               pdb_relative=anchor['pdb_relative'], pdb_sha256=anchor['pdb_sha256'], anchor_design=request['anchor_design'],
               added_since_anchor=';'.join(new_added), historical_scan_branchwidth=scan[delta]['branchwidth'] if delta>=0 else 'unmeasured',
               destination=request['destination'], kind=request['kind'])
    designs.append(row)
write(PACKAGE/'designs.tsv', designs)
(PACKAGE/'designs.json').write_text(json.dumps(designs, indent=2)+'\n')
paper = json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']
mark = {r['design_id']: r for r in paper['runs'] if r['method']=='MARK*'}
audits = []
for row in designs:
    name = row['design_id']
    target = PACKAGE/'preflight'/name
    target.mkdir()
    original = shlex.split((Path(mark[row['anchor_design']]['path'])/'command.sh').read_text())
    mark_cp = original[original.index('-cp')+1]
    assert all(Path(p).exists() for p in mark_cp.split(':'))
    command = [JAVA, '--add-opens', 'java.base/java.util=ALL-UNNAMED', '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
               '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED', '-Xmx8g', '-XX:ActiveProcessorCount=4',
               '-Djava.io.tmpdir='+os.environ['TMPDIR'], '-Dosprey.bench.method=sequence_dump', '-Dosprey.bench.numCPUs=4',
               '-Dosprey.bench.numGpus=0', '-Dosprey.bench.designId='+name, '-Dosprey.bench.outputDir='+str(target),
               '-Dosprey.bench.pdbPath='+str(PACKAGE/row['pdb_relative']), '-Dosprey.bench.mutable='+row['mutable'],
               '-Dosprey.bench.flexible='+row['flexible'], '-cp', mark_cp, 'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    (target/'command.sh').write_text(shlex.join(command)+'\n')
    with (target/'preflight.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    log = (target/'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log
    positions = re.findall(r'\[FRONTIER_POSITION\] state=Complex residue=(\S+)', log)
    assert len(positions)==row['total_positions']
    assert set(positions)==set((row['mutable']+';'+row['flexible']).split(';'))
    sequences = read(target/(name+'_sequences.tsv'))
    anchor_sequences = read(BASE/'preflight'/row['anchor_design']/(row['anchor_design']+'_sequences.tsv'))
    assert len(sequences)==row['expected_sequences']
    assert [r['sequence'] for r in sequences]==[r['sequence'] for r in anchor_sequences]
    env = dict(os.environ, JAVA_HOME=str(Path(JAVA).parent.parent), BUILD_ROOT=str(build), INPUT_ROOT=str(PACKAGE),
               RESULT_ROOT=str(ROOT/'pack_preflight'), PACKSTAR_EXPECT_GPU='RTX PRO 6000')
    with (target/'pack_preflight_driver.log').open('w') as log:
        subprocess.run(['python3', str(ROOT/'portable_run.py'), '--frontier-manifest', str(PACKAGE/'designs.tsv'),
                        '--index', str(row['task_id']), '--mode', 'preflight', '--arm', 'pair-only', '--rb', '1',
                        '--seed', '42', '--gpus', '4', '--gpu-gib', '85', '--heap-gib', '8', '--host-gib', '4'],
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    manifests = list((ROOT/'pack_preflight').glob(name+'_preflight_*/manifest.json'))
    assert len(manifests)==1
    result = json.loads(manifests[0].read_text())
    assert result['status']=='COMPLETED' and result['result_rows']==row['expected_sequences']
    pack_log = (manifests[0].parent/'run.log').read_text()
    assert int(re.search(r'Complex positions: (\d+)', pack_log)[1])==row['total_positions']
    audits.append(dict(design_id=name, positions=row['total_positions'], sequences=len(sequences),
                       anchor_sequences_identical=True, mark_preflight='PASS', pack_preflight='PASS'))
    print('PREFLIGHT_OK', json.dumps(audits[-1]), flush=True)

protocol = dict(designs=3, cpus=128, gpus=4, gpu_model='RTX PRO 6000', heap_gib=850, host_gib=800, gpu_gib=85,
                seed=42, rb=1, arm='pair-only', epsilon=0.683, global_limit_seconds=14400,
                per_design_max_seconds=3600, future_case_reserve_seconds=1200, shutdown_reserve_seconds=300,
                pack_under30_minutes='Evaluate after measurement; one-hour cap retains comparable timeout evidence.',
                source_build=str(VALIDATED/'build'), benchmark_mark_jobs_submitted=False)
(PACKAGE/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
(ROOT/'verification.json').write_text(json.dumps(dict(preflight=audits, protocol=protocol,
                 designs_sha256=digest(PACKAGE/'designs.tsv')), indent=2)+'\n')
(PACKAGE/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(PACKAGE))+'\n'
                                     for p in sorted(PACKAGE.rglob('*')) if p.is_file()))
write(HERE/'designs.tsv', designs)
registry = dict(prep_job=os.environ['SLURM_JOB_ID'], prep_root=str(ROOT), package=str(PACKAGE), designs=3,
                verification=audits, packstar_status='READY_TO_SUBMIT', markstar_status='PLAN_ONLY')
(HERE/'registry.json').write_text(json.dumps(registry, indent=2)+'\n')
(ROOT/'READY').write_text('Three designs passed MARK* and PACK* input preflights.\n')
print('READY', ROOT, flush=True)
