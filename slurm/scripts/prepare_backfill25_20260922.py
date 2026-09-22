"""Freeze exactly the 25 approved MARK* extensions and preflight them on Slurm."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

assert os.environ.get('SLURM_JOB_ID'), 'Submit through Slurm'
root = Path(sys.argv[1])
package = root / 'package'
repo = Path('/home/users/lz280/OSPREY3-fresh-packstar')
prior = Path('/usr/xtmp/lz280/packstar_flex_frontier38_20260920/prep_12651495')
frontier = Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config')
approved = {
    '2p4a': [7, 8, 9], '2rfd': [6, 7, 8], '4z80': [0, 1],
    '5em2': [5, 6], '2hnv': [4, 5], '2q2a': [2, 3],
    '1gwc': [2, 3], '2rl0': [1, 2], '2rfe': [2], '2hnu': [5],
    '2rf9': [11], '2xxm': [3], '3eb6': [6], '1b6c': [2], '3bu8': [11],
}
manifests = ['frontier.tsv', 'frontier_add4.tsv', 'frontier_small4.tsv',
             'frontier_remaining18_plus4.tsv']


def tsv(path):
    return list(csv.DictReader(io.StringIO(path.read_text().replace('\\t', '\t')), delimiter='\t'))


def dump(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


for name in ['structures', 'config', 'source', 'preflight']:
    (package / name).mkdir(parents=True, exist_ok=True)
assert (prior / 'READY').is_file()
cp = (prior / 'test_classpath.txt').read_text().strip()
(root / 'test_classpath.txt').write_text(cp + '\n')
shutil.copy2(prior / 'package/source/GenericPDBBench.java', package / 'source/GenericPDBBench.java')
assert '.setStabilityThreshold(null)' in (package / 'source/GenericPDBBench.java').read_text()
spec_path = Path('/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv')
specs = {r[0]: r for r in csv.reader(spec_path.open()) if r and not r[0].startswith('#')}
shutil.copy2(spec_path, package / 'config/base_specs.csv')
current = []
for name in manifests:
    shutil.copy2(repo / 'slurm/h200' / name, package / 'config' / name)
    current.extend(tsv(package / 'config' / name))
assert len(current) == 139 and len({r['design_id'] for r in current}) == 139
(package / 'config/approved_levels.json').write_text(json.dumps(approved, indent=2) + '\n')
pdbs = {s: Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped') / s / f'{s}.min.reduce.renum.pdb' for s in approved}
print('STORAGE_ESTIMATE', json.dumps(dict(pdb_files=len(pdbs), pdb_bytes=sum(p.stat().st_size for p in pdbs.values()),
      designs=25, preflight_directories=25, note='Input PDBs and small manifests only; search outputs remain in xtmp.')), flush=True)
designs = []
for system, levels in approved.items():
    spec = specs[system]
    mutable = spec[5].split(';')
    base = spec[6].split(';')
    if system == '2rfd':
        assert 'B552' in base
        base.remove('B552')  # Preserve the approved normalization from the previous batch.
    scan_path = frontier / 'expansion_scans' / f'{system}_scan.tsv'
    scans = {int(r['k']): r['added_residues'].split(';') if r['added_residues'] else [] for r in tsv(scan_path)}
    shutil.copy2(scan_path, package / 'config' / scan_path.name)
    removal = [r for r in reversed(base) if r[0] in {m[0] for m in mutable}]

    def flex_at(k):
        return [r for r in base if r not in removal[:-k]] if k < 0 else base + scans[k]

    siblings = [r for r in current if r['system'] == system]
    assert siblings and max(int(r['flex_delta']) for r in siblings) + 1 == min(levels), system
    for sibling in siblings:
        assert sibling['mutable'] == ';'.join(mutable), sibling
        assert sibling['flexible'] == ';'.join(flex_at(int(sibling['flex_delta']))), sibling
    pdb = package / 'structures' / pdbs[system].name
    shutil.copy2(pdbs[system], pdb)
    checksum = hashlib.sha256(pdb.read_bytes()).hexdigest()
    assert all(r['pdb_sha256'] == checksum for r in siblings), system
    for k in levels:
        flex = flex_at(k)
        assert len(flex) == len(base) + k
        assert len(set(mutable + flex)) == len(mutable + flex)
        did = f'{system}_flex_p{k}'
        assert did not in {r['design_id'] for r in current}
        designs.append(dict(task_id=len(designs), design_id=did, system=system, flex_delta=k,
            mutable=';'.join(mutable), flexible=';'.join(flex), added=';'.join(scans[k]), removed='',
            wt_flex_count=len(flex), total_positions=len(mutable + flex),
            expected_sequences=int(siblings[0]['expected_sequences']), pdb_relative='structures/' + pdb.name,
            pdb_sha256=checksum, markstar_cpus=16, markstar_memory_gib=192, markstar_heap_gib=160))
assert len(designs) == 25 and '3k3q' not in approved
java = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'
seqsets = {}
for row in designs:
    out = package / 'preflight' / row['design_id']
    out.mkdir()
    command = [java, '--add-opens', 'java.base/java.util=ALL-UNNAMED', '--add-opens',
        'java.base/java.lang=ALL-UNNAMED', '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED',
        '-Xmx8g', '-XX:ActiveProcessorCount=4', '-Djava.io.tmpdir=' + os.environ['TMPDIR'],
        '-Dosprey.bench.method=sequence_dump', '-Dosprey.bench.numCPUs=4', '-Dosprey.wmb.numGpus=0',
        '-Dosprey.bench.designId=' + row['design_id'], '-Dosprey.bench.outputDir=' + str(out),
        '-Dosprey.bench.pdbPath=' + str(package / row['pdb_relative']),
        '-Dosprey.bench.mutable=' + row['mutable'], '-Dosprey.bench.flexible=' + row['flexible'],
        '-cp', cp, 'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    with (out / 'preflight.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
    log = (out / 'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in log and 'WARNING: mutable residue' not in log, row
    positions = [x.split(' residue=')[1].split()[0] for x in log.splitlines() if x.startswith('[FRONTIER_POSITION] state=Complex ')]
    assert len(positions) == row['total_positions'] and set(positions) == set((row['mutable'] + ';' + row['flexible']).split(';')), row
    seqs = tsv(out / (row['design_id'] + '_sequences.tsv'))
    assert len(seqs) == row['expected_sequences'], row
    values = [r['sequence'] for r in seqs]
    system = row['system']
    if system in seqsets:
        assert seqsets[system] == values, row
    else:
        sibling = next(r for r in current if r['system'] == system)
        source = Path(sibling['source_package']) / 'preflight' / sibling['design_id'] / (sibling['design_id'] + '_sequences.tsv')
        assert values == [r['sequence'] for r in tsv(source)], row
        seqsets[system] = values
    print('PREFLIGHT_OK', row['task_id'], row['design_id'], row['total_positions'], len(seqs), flush=True)
(package / 'designs.json').write_text(json.dumps(designs, indent=2) + '\n')
dump(package / 'designs.tsv', designs)
protocol = json.loads((prior / 'package/protocol.json').read_text())
protocol.update(memory_gib=192, heap_gib=160, num_cpus=16, account='grisman', partition='compsci',
    exclude_nodes='', time_limit_days=14, resources_by_system={s: dict(memory_gib=192, heap_gib=160) for s in approved},
    selection='25 explicitly approved upward extensions on 15 systems; 3K3Q +1/+2 excluded.',
    level_choice='Exact approved levels; frozen expansion order verified against all existing sibling configurations.',
    status='25 actual position and sequence preflights passed; ready for production.',
    compiled_protocol_prep=str(prior), prediction='Memory allocation and previous runtimes do not guarantee completion.')
(package / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
(package / 'README.md').write_text('# MARK* 25-design extension\n\n16 CPU, 192 GiB memory, 160 GiB heap, compsci, account=grisman, 14 days.\n25 exact position and sequence preflights passed; all previous sibling configurations and sequence orders matched.\n2RFD continues to exclude B552. No 3K3Q extension.\n')
checks = [hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + str(p.relative_to(package)) for p in sorted(package.rglob('*')) if p.is_file()]
(package / 'SHA256SUMS').write_text('\n'.join(checks) + '\n')
(root / 'READY').write_text('25 designs validated\n')
print('READY', root, flush=True)
