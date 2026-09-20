"""Freeze the approved four-hour-cutoff addition and run actual sequence preflights."""
import csv, hashlib, io, json, os, shutil, subprocess, sys
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
root=Path(sys.argv[1]); package=root/'package'
for name in ['structures','preflight','config','source']:
    (package/name).mkdir(parents=True,exist_ok=True)
draft=Path('/usr/xtmp/lz280/frontier_review_12651281/draft_designs.tsv')
exclude={'4z80','2rl0','2rfe','2q2a','1gwc','1b6c'}
rows=[r for r in csv.DictReader(draft.open(),delimiter='\t') if not (r['tier']=='high' and r['system'] in exclude)]
assert len(rows)==52 and len({r['design_id'] for r in rows})==52
prior=Path('/usr/xtmp/lz280/packstar_flex_frontier20_20260919/prep_12643020')
assert (prior/'READY').is_file()
cp=(prior/'test_classpath.txt').read_text().strip()
(root/'test_classpath.txt').write_text(cp+'\n')
shutil.copy2(prior/'package/source/GenericPDBBench.java',package/'source/GenericPDBBench.java')
assert '.setStabilityThreshold(null)' in (package/'source/GenericPDBBench.java').read_text()
shutil.copy2(draft,package/'config/original_review.tsv')
systems=sorted({r['system'] for r in rows})
pdbs={s:Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped')/s/f'{s}.min.reduce.renum.pdb' for s in systems}
print('STORAGE_ESTIMATE',json.dumps(dict(pdb_files=len(pdbs),pdb_bytes=sum(p.stat().st_size for p in pdbs.values()),designs=52,preflight_directories=52,notes='Small input package only; production EMAT/logs remain in xtmp.')),flush=True)
for s,p in pdbs.items():
    shutil.copy2(p,package/'structures'/p.name)
    shutil.copy2(Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config/expansion_scans')/f'{s}_scan.tsv',package/'config'/f'{s}_scan.tsv')
designs=[]; seqsets={}; normalizations=[]; failures=[]
java='/home/users/lz280/java/jdk-17.0.2+8/bin/java'
groups={'compsci':[],'grisman':[],'fennario':[]}
for i,r in enumerate(rows):
    s=r['system']; mut=r['mutable'].split(';'); flex=r['flexible'].split(';')
    if s=='3cal':
        old_id=r['design_id']; d=int(r['flex_delta'])
        scans={int(x['k']):x for x in csv.DictReader(io.StringIO((package/'config/3cal_scan.tsv').read_text().replace('\\t','\t')),delimiter='\t')}
        extra=[x for x in scans[d+1]['added_residues'].split(';') if x not in scans[d]['added_residues'].split(';')]
        assert 'D210' in flex and len(extra)==1
        flex=[x for x in flex if x!='D210']+extra
        r.update(flexible=';'.join(flex),flex_delta=str(d+1),design_id=f'3cal_flex_p{d+1}')
        normalizations.append(dict(original_design_id=old_id,design_id=r['design_id'],removed_invalid='D210',added=extra[0],reason='D210 automatically deleted by OSPREY template matching; preserve approved actual position count using next frozen scan residue. Normalized baseline is 10 positions.'))
    assert len(set(mut+flex))==int(r['total_positions'])
    group='fennario' if r['tier'] in ['low','middle'] else 'compsci' if r['tier']=='extra_plus1' or s in ['2hnu','2hnv'] else 'grisman'
    groups[group].append(i)
    row=dict(task_id=i,design_id=r['design_id'],system=s,flex_delta=int(r['flex_delta']),mutable=r['mutable'],flexible=r['flexible'],wt_flex_count=len(flex),total_positions=len(mut+flex),expected_sequences=21 if s=='4wem' else 1+19*len(mut),pdb_relative='structures/'+pdbs[s].name,pdb_sha256=hashlib.sha256((package/'structures'/pdbs[s].name).read_bytes()).hexdigest(),markstar_memory_gib=96,markstar_heap_gib=64,tier=r['tier'],scheduling_group=group,review_index=int(r['review_index']))
    designs.append(row)
    out=package/'preflight'/row['design_id'];out.mkdir()
    cmd=[java,'--add-opens','java.base/java.util=ALL-UNNAMED','--add-opens','java.base/java.lang=ALL-UNNAMED','--add-opens','java.base/java.lang.invoke=ALL-UNNAMED','-Xmx8g','-XX:ActiveProcessorCount=4','-Djava.io.tmpdir='+os.environ['TMPDIR'],'-Dosprey.bench.method=sequence_dump','-Dosprey.bench.numCPUs=4','-Dosprey.wmb.numGpus=0','-Dosprey.bench.designId='+row['design_id'],'-Dosprey.bench.outputDir='+str(out),'-Dosprey.bench.pdbPath='+str(package/row['pdb_relative']),'-Dosprey.bench.mutable='+row['mutable'],'-Dosprey.bench.flexible='+row['flexible'],'-cp',cp,'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    with (out/'preflight.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=240)
    log=(out/'preflight.log').read_text()
    if 'WARNING: flexible residue' in log or 'WARNING: mutable residue' in log:
        failures.append(row); print('PREFLIGHT_INVALID_RESIDUE',row,flush=True); continue
    positions=[x for x in log.splitlines() if x.startswith('[FRONTIER_POSITION] state=Complex ')]
    assert len(positions)==row['total_positions'] and {x.split(' residue=')[1].split()[0] for x in positions}==set(mut+flex),row
    seqs=list(csv.DictReader((out/(row['design_id']+'_sequences.tsv')).open(),delimiter='\t'))
    assert len(seqs)==row['expected_sequences'],row
    vals=[v['sequence'] for v in seqs]
    if s in seqsets: assert seqsets[s]==vals,row
    seqsets[s]=vals
    print('PREFLIGHT_OK',i,row['design_id'],row['total_positions'],len(seqs),group,flush=True)
assert {k:len(v) for k,v in groups.items()}=={'compsci':6,'grisman':10,'fennario':36}
assert not failures,failures
(package/'config/input_normalizations.json').write_text(json.dumps(normalizations,indent=2)+'\n')
(package/'designs.json').write_text(json.dumps(designs,indent=2)+'\n')
with (package/'designs.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(designs[0]),delimiter='\t');w.writeheader();w.writerows(designs)
(root/'groups.json').write_text(json.dumps(groups,indent=2)+'\n')
protocol=json.loads((prior/'package/protocol.json').read_text())
protocol.update(resources_by_system={s:dict(memory_gib=96,heap_gib=64) for s in systems},partition='group-specific: compsci, grisman non-fennario, grisman fennario-only',exclude_nodes='group-specific',selection='User-approved four-hour cutoff: omit six original high tiers; retain all 36 low/middle, 12 high and four extra +1.',level_choice='Exact approved review TSV with six high-tier exclusions.',time_limit_days=14,compiled_protocol_prep=str(prior),status='52 preflighted designs; no production job yet',prediction='3x per position heuristic, not a completion guarantee; 64 GiB heap not certified by sequence preflight')
(package/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
(package/'README.md').write_text('# Approved 52-design addition\n\n52 actual position and sequence preflights passed. 16 CPU, 96 GiB memory, 64 GiB heap, FP64 CPU CCD, epsilon 0.683, stability filter disabled, fresh EMAT per design, 14-day limit. Separate batch; prior frozen artifacts unchanged.\n')
checks=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(package)) for p in sorted(package.rglob('*')) if p.is_file()]
(package/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(root/'READY').write_text('52 designs validated\n')
print('GROUPS',json.dumps(groups),flush=True)
print('READY',root,flush=True)
