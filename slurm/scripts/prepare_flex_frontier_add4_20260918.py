"""Freeze and validate four additions, without changing the running twelve-system cohort."""
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
SYSTEMS = {'2rfd': list(range(5)), '3u7y': list(range(7)), '3gxu': list(range(10)), '3bua': [-2,-1,0]}

PREVIOUS_PREP = Path('/usr/xtmp/lz280/packstar_flex_frontier10_20260917/prep_12631844')

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
          'notes':'Four PDBs, small manifests and one Java overlay. No dataset download. Long-run caches/log sizes depend on search; all stay in xtmp.'}
(ROOT/'storage_estimate.json').write_text(json.dumps(estimate,indent=2)+'\n')
print('INPUT_ESTIMATE',json.dumps(estimate),flush=True)

designs=[]
normalizations=[]
for system,levels in SYSTEMS.items():
    spec=specs[system]; mutable=spec[5].split(';'); rawflex=spec[6].split(';')
    baseflex=list(dict.fromkeys(rawflex))
    if system == '2rfd':
        assert 'B552' in baseflex
        baseflex.remove('B552')
        normalizations.append(dict(system=system, original_flexible=rawflex,
            normalized_flexible=baseflex,
            reason='B552 is automatically deleted by OSPREY template matching in the frozen PDB; exclude the inactive requested site explicitly. Effective baseline is 8 positions, not the historical list count 9.'))
    if rawflex!=baseflex and system != '2rfd':
        assert system=='4znc' and rawflex.count('F731')==2
        normalizations.append(dict(system=system,original_flexible=rawflex,normalized_flexible=baseflex,
                                   reason='Remove duplicate residue F731; retain original first-occurrence order'))
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
        if (14.0 if system == '3bua' else float(screen[system]['mark_days']))>=7: assert k<=0
        did=f'{system}_flex_{"m"+str(-k) if k<0 else "p"+str(k)}'
        row=dict(task_id=len(designs),design_id=did,system=system,flex_delta=k,
                 mutable=';'.join(mutable),flexible=';'.join(flex),
                 added=';'.join(added),removed=';'.join(removed),
                 wt_flex_count=len(flex),total_positions=len(flex)+len(mutable),
                 expected_sequences=1+19*len(mutable),
                 pdb_relative='structures/'+pdb.name,pdb_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
                 historical_mark_days=(14.0 if system == '3bua' else float(screen[system]['mark_days'])),
                 historical_pack_minutes=float(screen[system]['sep05_seconds'])/60 if screen[system]['sep05_seconds'] else None,
                 historical_mark_kind='timeout_lower_bound' if system == '3bua' else 'result_csv',
                 markstar_memory_gib=96 if system == '3bua' else 128,
                 markstar_heap_gib=64 if system == '3bua' else 96,
                 rough_mark_days_2x=(14.0 if system == '3bua' else float(screen[system]['mark_days']))*2**k,
                 rough_mark_days_3x=(14.0 if system == '3bua' else float(screen[system]['mark_days']))*3**k,
                 old_scan_branchwidth=scans[k]['branchwidth'] if k>=0 else 'unmeasured')
        designs.append(row)
assert len(designs)==25 and len(SYSTEMS)==4
dump_tsv(package/'designs.tsv',designs)
(package/'designs.json').write_text(json.dumps(designs,indent=2)+'\n')
(package/'config/input_normalizations.json').write_text(json.dumps(normalizations,indent=2)+'\n')
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
assert source==(PREVIOUS_PREP/'package/source/GenericPDBBench.java').read_text(), 'MARK* protocol overlay changed'
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
    'resources_by_system': {s:dict(memory_gib=96 if s=='3bua' else 128, heap_gib=64 if s=='3bua' else 96) for s in SYSTEMS},
    'time_limit_days':14,'account':'grisman','partition':'compsci',
    'markstar_reduce_minimizations':True,'markstar_correction_tightening':True,
    'markstar_full_parallel_leaf_batch':False,
    'energy_matrices':'fresh per design, no historical or sibling cache reuse',
    'selection':'Four user-selected systems supplement the existing 12 systems / 46 designs; no existing design is replaced',
    'level_choice':'2rfd 0..4, 3u7y 0..6, 3gxu 0..9: nearest integer log_3(14 days / historical MARK time), including every intermediate level. These are heuristic search boundaries, not measured frontiers. 3bua -2,-1,0 explicitly selected by user.',
    'removal_rule':'Reverse original protein-side WT-flex order, keep all ligand-side flexibility and all mutable sites',
    'historical_3bua':'MARK job 11899111 TIMEOUT at 14 days, MaxRSS 39681640 KiB, Xmx64g; not OOM. PACK A12509108 T15 failed host budget check, a separate cause.',
    'prediction':'2x-3x growth per site is heuristic; disabling historical stability filtering and hardware differences affect runtime',
    'timeouts':'Right-censored at 14 days, not completed timing',
    'packstar_status':'Additional H200 cohort only; preserve the running 12-system build, inputs, indices and outputs',
}
(package/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
lines=['# Four additions: 25 designs, combined cohort 16 systems / 71 designs','',
       '| System | Flex delta | Total positions | Memory GiB | Heap GiB |',
       '|---|---|---|---|---|']
for system, levels in SYSTEMS.items():
    group=[r for r in designs if r['system']==system]
    lines.append(f"| {system} | {levels} | {[r['total_positions'] for r in group]} | {group[0]['markstar_memory_gib']} | {group[0]['markstar_heap_gib']} |")
lines += ['', 'All 25 actual OSPREY position/sequence preflights passed. Fixed mutable residues and sequence order per system.',
          '16 CPU, FP64, CPU CCD, epsilon=0.683, stability filter off, fresh EMAT, 14-day limit.',
          'Original 12-system H200 experiments stay intact. Use a separate additional input package and result root.']
(package/'README.md').write_text('\n'.join(lines)+'\n')
checks=[]
for path in sorted(package.rglob('*')):
    if path.is_file(): checks.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(package)))
(package/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(ROOT/'READY').write_text('25 designs validated\n')
print('\n'.join(lines),flush=True)
print('PREP_ROOT='+str(ROOT),flush=True)
