"""Freeze eight long-run targets and one conservative tier with matching inputs."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Submit preparation through Slurm'
HERE = Path(__file__).parent
if (HERE/'registry.json').exists():
    assert not json.loads((HERE/'registry.json').read_text()).get('submissions'), 'Already submitted; do not replay preparation'
REPO = HERE.parents[2]
ROOT = Path('/usr/xtmp/lz280/markstar_near14_20261006') / ('prep_' + os.environ['SLURM_JOB_ID'])
ROOT.mkdir(parents=True, exist_ok=False)
PACKAGE = ROOT / 'package'
BASE_INPUT = Path('/usr/xtmp/lz280/markstar_backfill_20260922/handoff_12686177/package')
VALIDATED_RUN = Path('/usr/xtmp/lz280/packstar_gpu101_20261005/run_12814879')
CPU_RUN = Path('/usr/xtmp/lz280/packstar_cpu32_optimized_20261004/run_12775476')
SCAN_ROOT = Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config/expansion_scans')
JAVA = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'

def read_tsv(path):
    return list(csv.DictReader(io.StringIO(path.read_text().replace('\\t', '\t')), delimiter='\t'))

def write_tsv(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

requested = read_tsv(HERE/'requested_designs.tsv')
gpu_cases = {r['design'] for r in read_tsv(VALIDATED_RUN/'cases.tsv')}
cpu_cases = {r['design'] for r in read_tsv(CPU_RUN/'cases.tsv')}
assert (CPU_RUN/'run.py').read_bytes() == (VALIDATED_RUN/'run.py').read_bytes()
assert (CPU_RUN/'build/test_classpath.txt').read_bytes() == (VALIDATED_RUN/'build/test_classpath.txt').read_bytes()
assert (CPU_RUN/'build/source/slurm/h200/production.properties').read_bytes() == (VALIDATED_RUN/'build/source/slurm/h200/production.properties').read_bytes()
compile((HERE/'run_group.py').read_text(), str(HERE/'run_group.py'), 'exec')
subprocess.run(['bash', '-n', str(HERE/'run_group.slurm')], check=True)
assert requested and len({r['design_id'] for r in requested}) == len(requested)
prior_rows = {r['design_id']: r for r in read_tsv(BASE_INPUT/'designs.tsv')}
assert not set(prior_rows) & {r['design_id'] for r in requested}
paper = json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']
mark = {r['design_id']: r for r in paper['runs'] if r['method']=='MARK*'}
structures = {r['system']: BASE_INPUT/prior_rows[r['anchor_design']]['pdb_relative'] for r in requested}
estimate = dict(input_pdb_files=len(structures), input_pdb_bytes=sum(p.stat().st_size for p in structures.values()),
                preparation_expected_files_under=150, preparation_expected_bytes_under=50*2**20,
                production_expected=f'{len(requested)} workloads; ordinary logs/caches expected in tens of GiB, with inherited DP table budgets. No dataset or build copy.',
                inherited_per_design_mapped_table_budget_gib=512, storage_root=str(ROOT.parent))
print('STORAGE_ESTIMATE', json.dumps(estimate), flush=True)
(ROOT/'storage_estimate.json').write_text(json.dumps(estimate, indent=2)+'\n')
for name in ['structures', 'preflight', 'config', 'provenance']:
    (PACKAGE/name).mkdir(parents=True)
for system, source in structures.items():
    shutil.copy2(source, PACKAGE/'structures'/source.name)
    shutil.copy2(SCAN_ROOT/(system+'_scan.tsv'), PACKAGE/'config'/(system+'_scan.tsv'))
shutil.copy2(HERE/'requested_designs.tsv', PACKAGE/'provenance/requested_designs.tsv')
shutil.copy2(__file__, ROOT/'prepare.py')
shutil.copy2(HERE/'run_group.py', ROOT/'run_group.py')
shutil.copy2(HERE/'run_group.slurm', ROOT/'run_group.slurm')

# Reuse the exact validated optimized build used by the current PRO6000 run.
old_build = VALIDATED_RUN/'build'
cp = (old_build/'test_classpath.txt').read_text().strip()
assert all(Path(p).exists() for p in cp.split(':'))
build = ROOT/'build'
(build/'source/slurm/h200').mkdir(parents=True)
for name in ['test_classpath.txt', 'git_head.txt', 'source.sha256']:
    shutil.copy2(old_build/name, build/name)
shutil.copy2(old_build/'source/slurm/h200/production.properties', build/'source/slurm/h200/production.properties')
shutil.copy2(VALIDATED_RUN/'run.py', ROOT/'portable_run.py')
candidate = Path(json.loads((VALIDATED_RUN/'protocol.json').read_text())['optimized_build'])
with (ROOT/'source_check.log').open('w') as f:
    subprocess.run(['diff', '-qr', str(REPO/'src/main'), str(candidate/'source/src/main')], stdout=f, stderr=subprocess.STDOUT, check=True)

designs = []
for index, request in enumerate(requested):
    anchor = prior_rows[request['anchor_design']]
    system = request['system']
    gpu_anchors = sorted(n for n in gpu_cases if n.startswith(system+'_flex_'))
    cpu_anchors = sorted(n for n in cpu_cases if n.startswith(system+'_flex_'))
    assert gpu_anchors or cpu_anchors, system
    # User explicitly moved the two CPU-cohort extensions to GPU as well.
    platform = 'gpu'
    reference_platform = 'gpu' if gpu_anchors else 'cpu'
    assert request['anchor_design'] in (gpu_cases if gpu_anchors else cpu_cases)
    tier, previous_tier = int(request['flex_delta']), int(anchor['flex_delta'])
    scans = {int(r['k']): r for r in read_tsv(PACKAGE/'config'/(system+'_scan.tsv'))}
    mutable = anchor['mutable'].split(';')
    if tier < 0:
        assert system == '3bua' and tier == -4 and previous_tier == -2
        baseline = prior_rows['3bua_flex_p0']['flexible'].split(';')
        protein_chains = {r[0] for r in mutable}
        removal = [r for r in reversed(baseline) if r[0] in protein_chains]
        assert anchor['removed'].split(';') == removal[:2]
        removed = removal[:4]
        flexible = [r for r in baseline if r not in removed]
        assert set(flexible) < set(anchor['flexible'].split(';'))
        all_added, new_added = [], []
    else:
        old_added = [x for x in scans[previous_tier]['added_residues'].split(';') if x]
        all_added = [x for x in scans[tier]['added_residues'].split(';') if x]
        assert all_added[:len(old_added)] == old_added
        new_added = all_added[len(old_added):]
        assert len(new_added) == tier-previous_tier
        flexible = anchor['flexible'].split(';') + new_added
        removed = [r for r in anchor.get('removed', '').split(';') if r]
    assert len(set(mutable+flexible)) == len(mutable+flexible)
    assert len(mutable+flexible) == int(anchor['total_positions']) + tier-previous_tier
    pdb = PACKAGE/anchor['pdb_relative']
    assert digest(pdb) == anchor['pdb_sha256']
    pdb_residues = {line[21].strip()+line[22:26].strip()+line[26].strip()
                    for line in pdb.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))}
    assert set(mutable+flexible) <= pdb_residues
    row = dict(task_id=index, design_id=request['design_id'], system=system, flex_delta=tier,
               mutable=anchor['mutable'], flexible=';'.join(flexible), added=';'.join(all_added),
               removed=';'.join(removed), wt_flex_count=len(flexible), total_positions=len(mutable+flexible),
               expected_sequences=int(anchor['expected_sequences']), pdb_relative=anchor['pdb_relative'], pdb_sha256=anchor['pdb_sha256'],
               anchor_design=request['anchor_design'], added_since_anchor=';'.join(new_added),
               historical_scan_branchwidth=scans[tier]['branchwidth'] if tier >= 0 else 'unmeasured',
               destination=request['destination'], kind=request['kind'], pack_platform=platform,
               pack_reference_job='12814879' if reference_platform == 'gpu' else '12775476',
               pack_reference_platform=reference_platform)
    designs.append(row)
write_tsv(PACKAGE/'designs.tsv', designs)
(PACKAGE/'designs.json').write_text(json.dumps(designs, indent=2)+'\n')

audits = []
for row in designs:
    name = row['design_id']
    target = PACKAGE/'preflight'/name
    target.mkdir()
    old = Path(mark[row['anchor_design']]['path'])
    old_command = shlex.split((old/'command.sh').read_text())
    mark_cp = old_command[old_command.index('-cp')+1]
    assert all(Path(p).exists() for p in mark_cp.split(':'))
    command = [JAVA, '--add-opens', 'java.base/java.util=ALL-UNNAMED', '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
               '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED', '-Xmx8g', '-XX:ActiveProcessorCount=4',
               '-Djava.io.tmpdir='+os.environ['TMPDIR'], '-Dosprey.bench.method=sequence_dump', '-Dosprey.bench.numCPUs=4',
               '-Dosprey.bench.numGpus=0', '-Dosprey.bench.designId='+name, '-Dosprey.bench.outputDir='+str(target),
               '-Dosprey.bench.pdbPath='+str(PACKAGE/row['pdb_relative']), '-Dosprey.bench.mutable='+row['mutable'],
               '-Dosprey.bench.flexible='+row['flexible'], '-cp', mark_cp, 'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    (target/'command.sh').write_text(shlex.join(command)+'\n')
    with (target/'preflight.log').open('w') as f:
        subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=240)
    log = (target/'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log, name
    actual_positions = re.findall(r'\[FRONTIER_POSITION\] state=Complex residue=(\S+)', log)
    assert len(actual_positions) == row['total_positions'], (name, actual_positions)
    assert set(actual_positions) == set((row['mutable']+';'+row['flexible']).split(';'))
    sequences = read_tsv(target/(name+'_sequences.tsv'))
    original_sequences = read_tsv(BASE_INPUT/'preflight'/row['anchor_design']/(row['anchor_design']+'_sequences.tsv'))
    assert len(sequences) == row['expected_sequences']
    assert [r['sequence'] for r in sequences] == [r['sequence'] for r in original_sequences]
    env = dict(os.environ, JAVA_HOME=str(Path(JAVA).parent.parent), BUILD_ROOT=str(build), INPUT_ROOT=str(PACKAGE),
               RESULT_ROOT=str(ROOT/'pack_preflight'), PACKSTAR_EXPECT_GPU='RTX PRO 6000')
    with (target/'pack_preflight_driver.log').open('w') as f:
        subprocess.run(['python3', str(ROOT/'portable_run.py'), '--frontier-manifest', str(PACKAGE/'designs.tsv'),
                        '--index', str(row['task_id']), '--mode', 'preflight', '--arm', 'pair-only', '--rb', '1',
                        '--seed', '42', '--gpus', '4' if row['pack_platform']=='gpu' else '0', '--gpu-gib', '85', '--heap-gib', '8', '--host-gib', '4'],
                       env=env, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=240)
    result_paths = list((ROOT/'pack_preflight').glob(name+'_preflight_*/manifest.json'))
    assert len(result_paths)==1
    result = json.loads(result_paths[0].read_text())
    assert result['status']=='COMPLETED' and result['result_rows']==row['expected_sequences']
    pack_log = (result_paths[0].parent/'run.log').read_text()
    assert int(re.search(r'Complex positions: (\d+)', pack_log)[1]) == row['total_positions']
    audits.append(dict(design_id=name, positions=row['total_positions'], sequences=len(sequences),
                       anchor_sequences_identical=True, mark_preflight='PASS', pack_preflight='PASS',
                       pack_platform=row['pack_platform'], pack_reference_job=row['pack_reference_job']))
    print('PREFLIGHT_OK', json.dumps(audits[-1]), flush=True)

protocol = dict(designs=len(designs), systems=len(structures), existing_designs_modified=False, source_build=str(candidate),
                inherited_packstar_jobs=[12814879,12775476], heap_gib=850,
                host_gib=800, gpu_gib=85, seed=42, rb=1, arm='pair-only', epsilon=0.683,
                scheduling='Nine GPU workloads in one allocation after the existing batch; append only',
                platforms={'gpu': dict(cpus=128,gpus=4,model='RTX PRO 6000',global_limit_seconds=21600)},
                per_design_max_seconds=3600,
                future_case_reserve_seconds=1200, shutdown_reserve_seconds=300,
                markstar_status='MARK* may start after input preflight, independently of PACK* measurements')
(PACKAGE/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
(ROOT/'verification.json').write_text(json.dumps(dict(preflight=audits, protocol=protocol, designs_sha256=digest(PACKAGE/'designs.tsv')), indent=2)+'\n')
checks = [digest(p)+'  '+str(p.relative_to(PACKAGE)) for p in sorted(PACKAGE.rglob('*')) if p.is_file()]
(PACKAGE/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
# Small reviewable registry; large artifacts remain under xtmp.
write_tsv(HERE/'new_designs.tsv', [r for r in designs if r['kind'] != 'retained_extension'])
shutil.copy2(PACKAGE/'designs.tsv', HERE/'pack_designs.tsv')
registry = dict(prep_job=os.environ['SLURM_JOB_ID'], prep_root=str(ROOT), package=str(PACKAGE),
                design_count=len(designs), new_design_count=sum(r['kind'] != 'retained_extension' for r in designs), designs_sha256=digest(PACKAGE/'designs.tsv'), verification=audits,
                packstar_status='PREFLIGHT_PASSED_AWAITING_SUBMISSION', markstar_status='PLAN_ONLY',
                mode='APPEND_ONLY_KEEP_SMALLER_DESIGNS_FIRST', earlier_unsubmitted_preps=['12827301','12827308'])
(HERE/'registry.json').write_text(json.dumps(registry, indent=2)+'\n')
(ROOT/'READY').write_text(f'{len(designs)} designs: both benchmark preflights passed\n')
print('READY', ROOT, flush=True)
