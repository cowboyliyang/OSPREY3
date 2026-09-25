"""Freeze and validate 32 nested designs on Slurm before MARK* submission."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

assert os.environ.get('SLURM_JOB_ID'), 'Run via Slurm'
ROOT = Path(sys.argv[1])
BUILD = Path('/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626762')
BASE = Path('/usr/xtmp/lz280/bench_comparison')
SCREEN = Path('/usr/xtmp/lz280/packstar_frontier_screen_20260917/12620972')
FRONTIER = Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2')
JAVA = '/home/users/lz280/java/jdk-17.0.2+8/bin/java'
JAVAC = '/home/users/lz280/java/jdk-17.0.2+8/bin/javac'
SYSTEMS = {
    '2xgy': [-2,-1,0], '4u3s': [-2,-1,0],
    '3ma2': [0,1], '4wyu': [0,1], '5a6y': [0,1,2],
    '3k3q': [0,1,2], '4z80': [0,1,2],
    '2rl0': [0,1,2,3], '2q2a': [0,1,2,3,4], '4wyq': [0,1,2,3],
}

def tsv(path):
    return list(csv.DictReader(io.StringIO(path.read_text().replace('\\t','\t')),delimiter='\t'))

def dump_tsv(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader(); w.writerows(rows)

ROOT.mkdir(parents=True,exist_ok=True)
package=ROOT/'package'
for name in ['structures','config','source','preflight']:
    (package/name).mkdir(parents=True,exist_ok=True)
(ROOT/'classes').mkdir(exist_ok=True)
specs={r[0]:r for r in csv.reader(BASE.joinpath('design_specs_prepped.csv').open()) if r and not r[0].startswith('#')}
screen={r['design']:r for r in tsv(SCREEN/'all38.tsv')}
inputs=[Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped')/s/f'{s}.min.reduce.renum.pdb' for s in SYSTEMS]
estimate={'pdb_files':len(inputs),'pdb_bytes':sum(p.stat().st_size for p in inputs),
          'designs':sum(map(len,SYSTEMS.values())),
          'notes':'Ten PDBs, small manifests and one Java overlay. No dataset download. Long-run caches/log sizes depend on search; all stay in xtmp.'}
(ROOT/'storage_estimate.json').write_text(json.dumps(estimate,indent=2)+'\n')
print('INPUT_ESTIMATE',json.dumps(estimate),flush=True)

designs=[]
for system,levels in SYSTEMS.items():
    spec=specs[system]; mutable=spec[5].split(';'); baseflex=spec[6].split(';')
    assert len(baseflex)==len(set(baseflex)) and not set(mutable)&set(baseflex)
    pdb=Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped')/system/f'{system}.min.reduce.renum.pdb'
    destination=package/'structures'/pdb.name
    shutil.copy2(pdb,destination)
    residues={}
    for line in pdb.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            resid=line[21].strip()+line[22:26].strip()+line[26].strip()
            residues.setdefault(resid,line[17:20].strip())
    scans={int(r['k']):r for r in tsv(FRONTIER/'config/expansion_scans'/f'{system}_scan.tsv')}
    shutil.copy2(FRONTIER/'config/expansion_scans'/f'{system}_scan.tsv',package/'config'/f'{system}_scan.tsv')
    # For the already-long systems, shrink protein-side WT flexibility in
    # reverse original-list order, retaining all ligand-side flexible sites.
    protein_chains={r[0] for r in mutable}
    removal=[r for r in reversed(baseflex) if r[0] in protein_chains]
    previous=set()
    for k in levels:
        removed=removal[:-k] if k<0 else []
        added=scans[k]['added_residues'].split(';') if k>0 else []
        flex=[r for r in baseflex if r not in removed]+added
        assert len(flex)==len(baseflex)+k and len(flex)==len(set(flex))
        assert not set(mutable)&set(flex) and set(mutable+flex)<=set(residues)
        assert previous<=set(flex); previous=set(flex)
        if float(screen[system]['mark_days'])>=7: assert k<=0
        did=f'{system}_flex_{"m"+str(-k) if k<0 else "p"+str(k)}'
        row=dict(task_id=len(designs),design_id=did,system=system,flex_delta=k,
                 mutable=';'.join(mutable),flexible=';'.join(flex),
                 added=';'.join(added),removed=';'.join(removed),
                 wt_flex_count=len(flex),total_positions=len(flex)+len(mutable),
                 expected_sequences=int(screen[system]['sep05_rows']),
                 pdb_relative='structures/'+pdb.name,pdb_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
                 historical_mark_days=float(screen[system]['mark_days']),
                 historical_pack_minutes=float(screen[system]['sep05_seconds'])/60,
                 rough_mark_days_2x=float(screen[system]['mark_days'])*2**k,
                 rough_mark_days_3x=float(screen[system]['mark_days'])*3**k,
                 old_scan_branchwidth=scans[k]['branchwidth'] if k>=0 else 'unmeasured')
        designs.append(row)
assert len(designs)==32 and len(SYSTEMS)==10
dump_tsv(package/'designs.tsv',designs)
(package/'designs.json').write_text(json.dumps(designs,indent=2)+'\n')
shutil.copy2(SCREEN/'all38.tsv',package/'config/historical_screening_all38.tsv')

# Small isolated protocol overlay: freeze the previously validated Java
# implementation, explicitly align stability filtering, and retain scores
# incrementally so partial results survive a two-week timeout.
original=BUILD/'source/src/test/java/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.java'
source=original.read_text()
needle='MARKStar.Settings.Builder sb = new MARKStar.Settings.Builder()\n                .setEpsilon(epsilon)'
assert source.count(needle)==1
source=source.replace(needle,needle+'\n                .setStabilityThreshold(null)\n                .addScoreFileWriter(new File(outputDir, designId + "_partial_scores.tsv"))')
needle='            case "sequence_dump":\n                runSequenceDump(confSpaces, outputDir, designId);'
assert source.count(needle)==1
source=source.replace(needle,needle+'''
                for (SimpleConfSpace space : new SimpleConfSpace[]{confSpaces.protein, confSpaces.ligand, confSpaces.complex}) {
                    String name = space == confSpaces.protein ? "Protein" : space == confSpaces.ligand ? "Ligand" : "Complex";
                    for (SimpleConfSpace.Position pos : space.positions) {
                        System.out.println("[FRONTIER_POSITION] state=" + name + " residue=" + pos.resNum + " rcs=" + pos.resConfs.size());
                    }
                }''')
overlay=package/'source/GenericPDBBench.java'
overlay.write_text(source)
shutil.copy2(original,package/'source/GenericPDBBench.original.java.txt')
cp=BUILD.joinpath('test_classpath.txt').read_text().strip()
subprocess.run([JAVAC,'-J-Xmx8g','-cp',cp,'-d',str(ROOT/'classes'),str(overlay)],check=True)
cp=str(ROOT/'classes')+':'+cp
(ROOT/'test_classpath.txt').write_text(cp+'\n')

sequence_sets={}
for row in designs:
    out=package/'preflight'/row['design_id']; out.mkdir()
    cmd=[JAVA,'--add-opens','java.base/java.util=ALL-UNNAMED','--add-opens','java.base/java.lang=ALL-UNNAMED',
         '--add-opens','java.base/java.lang.invoke=ALL-UNNAMED','-Xmx8g','-XX:ActiveProcessorCount=4',
         '-Djava.io.tmpdir='+os.environ['TMPDIR'],'-Dosprey.bench.method=sequence_dump',
         '-Dosprey.bench.numCPUs=4','-Dosprey.wmb.numGpus=0',
         '-Dosprey.bench.designId='+row['design_id'],'-Dosprey.bench.outputDir='+str(out),
         '-Dosprey.bench.pdbPath='+str(package/row['pdb_relative']),
         '-Dosprey.bench.mutable='+row['mutable'],'-Dosprey.bench.flexible='+row['flexible'],
         '-cp',cp,'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    with (out/'preflight.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
    text=(out/'preflight.log').read_text()
    assert 'WARNING: flexible residue' not in text and 'WARNING: mutable residue' not in text
    positions=[x for x in text.splitlines() if x.startswith('[FRONTIER_POSITION] state=Complex ')]
    assert len(positions)==row['total_positions'],row
    actual={x.split(' residue=')[1].split()[0] for x in positions}
    assert actual==set((row['mutable']+';'+row['flexible']).split(';')),row
    seqs=tsv(out/(row['design_id']+'_sequences.tsv'))
    assert len(seqs)==row['expected_sequences'],row
    sequence_values=[s['sequence'] for s in seqs]
    if row['system'] in sequence_sets: assert sequence_sets[row['system']]==sequence_values
    else: sequence_sets[row['system']]=sequence_values
    print('PREFLIGHT_OK',row['design_id'],'positions',len(positions),'sequences',len(seqs),flush=True)

protocol={
    'build_root':str(BUILD),'target_epsilon':0.683,'max_simultaneous_mutations':1,
    'stability_filter':False,'precision':'FP64','ccd':'CPU','num_cpus':16,
    'memory_gib':192,'heap_gib':160,'time_limit_days':14,'account':'grisman','partition_preference':'compsci then non-A5000 grisman',
    'markstar_reduce_minimizations':True,'markstar_correction_tightening':True,
    'markstar_full_parallel_leaf_batch':False,
    'energy_matrices':'fresh per design, no historical or sibling cache reuse',
    'timing':'retain end-to-end wall and separately reported search time; compare equivalent phases with PACK*',
    'selection':'10 historically favorable PACK* systems; this is a selected frontier cohort, not an unbiased 38-system aggregate',
    'prediction':'2x-3x per added WT flex is only a heuristic; hardware and removal of historical stability filtering may change runtime',
    'timeouts':'unfinished at 14 days are right-censored; do not report them as exact completed times',
    'packstar_status':'not submitted; use identical designs/PDB checksums/sequence manifests on future H200 runs',
}
(package/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
lines=['# Ten-system MARK* flexibility frontier','',
       '32 independent designs. No positive flex expansion for systems whose historical MARK* baseline exceeds seven days.',
       'Negative deltas remove protein-side WT flexible residues in reverse original-list order; ligand flexibility stays fixed.',
       'Positive deltas use the existing frozen expansion order. Mutable sites and enumerated sequence order stay fixed within each system.',
       'MARK* stability filtering is explicitly disabled to match PACK*. Historical runtimes are screening estimates, not predictions for the new protocol.','',
       '| System | Historical MARK* days | Flex deltas | Last-design rough days (2x–3x/site) |',
       '|---|---:|---|---:|']
for s,levels in SYSTEMS.items():
    rows=[r for r in designs if r['system']==s]; last=rows[-1]
    lines.append(f"| {s} | {last['historical_mark_days']:.2f} | {', '.join(f'{k:+d}' for k in levels)} | {last['rough_mark_days_2x']:.1f}–{last['rough_mark_days_3x']:.1f} |")
lines += ['', 'Execution: account=grisman, prefer partition=compsci, overflow to non-A5000 grisman; 16 CPUs, 192 GiB memory, 160 GiB Java heap, 14 days, no GPU.',
          'All designs are independently eligible for scheduling; no per-system or small-to-large dependencies.',
          'Input paths in designs.tsv are relative to this package so it can be copied to an H200 cluster. Build/classes/dependency paths are local and must be reconstructed there.',
          'All 32 OSPREY conformation-space preflights passed; positions and sequence counts/order were checked.',
          'No PACK* jobs submitted by this protocol.']
(package/'README.md').write_text('\n'.join(lines)+'\n')
checks=[]
for p in sorted(package.rglob('*')):
    if p.is_file(): checks.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(package)))
(package/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(ROOT/'READY').write_text('32 designs validated\n')
print('\n'.join(lines),flush=True)
print('PREP_ROOT='+str(ROOT),flush=True)
