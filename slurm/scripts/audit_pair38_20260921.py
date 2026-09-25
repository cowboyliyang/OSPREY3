"""Slurm-only baseline/resource census for the authorized full pair-only rerun."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get('SLURM_JOB_ID')
root = Path('/usr/xtmp/lz280/packstar_pair38_20260921') / ('audit_' + os.environ['SLURM_JOB_ID'])
root.mkdir(parents=True, exist_ok=False)
baseline = Path('/usr/xtmp/lz280/packstar_rb38_20260918/launch_12637884/rb1')
old_build = Path('/usr/xtmp/lz280/packstar_proposal_root_20260917/build_12632056')
new_build = Path('/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588')
rows = []
for path in sorted(baseline.glob('*/manifest.json')):
    m = json.loads(path.read_text())
    run = path.parent
    design = m['design']['design_id']
    result = run / (design + '_packstar.csv')
    records = list(csv.DictReader(result.open())) if result.exists() else []
    wall_text = (run / 'wall.time').read_text() if (run / 'wall.time').exists() else ''
    match = re.search(r'elapsed=([\d.]+) maxRssKiB=(\d+)', wall_text)
    p = m['properties']
    assert p['branchdp.cutoff.residualBudget'] == '1'
    assert p['packstar.pac.frequencySeverity.tripleEta'] == 'true'
    assert p['packstar.pac.frequencySeverity.jointMomentLearning'] == 'true'
    assert m['options']['seed'] == 42
    rows.append(dict(index=m['options']['index'], design=design,
        elapsed_s=float(match[1]) if match else m.get('end_epoch',m['start_epoch'])-m['start_epoch'],
        peak_rss_gib=int(match[2])/2**20 if match else None,
        status=m['status'], exit_code=m.get('exit_code'),
        rows=len(records), expected=m['expected_rows'], estimated=m.get('estimated_rows'),
        node=m['node'], gpu=(run/'gpus.csv').read_text().strip(), path=str(run)))
assert len(rows) == 38 and len({r['index'] for r in rows}) == 38
rows.sort(key=lambda r:r['elapsed_s'], reverse=True)
with (root / 'baseline.tsv').open('w') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(root/'baseline.json').write_text(json.dumps(rows,indent=2)+'\n')

diffs=[]
for subtree in ['src/main', 'slurm/h200/production.properties', 'slurm/h200/baseline38.csv']:
    a,b=old_build/'source'/subtree,new_build/'source'/subtree
    result=subprocess.run(['diff','-qr',str(a),str(b)],capture_output=True,text=True)
    diffs.append(dict(subtree=subtree,code=result.returncode,text=result.stdout+result.stderr))
(root/'source_comparison.json').write_text(json.dumps(diffs,indent=2)+'\n')
nodes=subprocess.check_output(['scontrol','show','nodes','--oneliner'],text=True)
(root/'scontrol_nodes.txt').write_text(nodes)
available=[]
for line in nodes.splitlines():
    d=dict(re.findall(r'(\w+)=([^ ]*)',line))
    partitions=set(d.get('Partitions','').split(','))
    if not partitions.intersection({'grisman','compsci-gpu'}): continue
    if any(s in d.get('State','') for s in ['DRAIN','DOWN','NOT_RESPONDING']): continue
    def tres(s):
        return dict(part.split('=',1) for part in s.split(',') if '=' in part)
    cfg,alloc=tres(d.get('CfgTRES','')),tres(d.get('AllocTRES',''))
    free_types={key[9:]:int(value)-int(alloc.get(key,0)) for key,value in cfg.items() if key.startswith('gres/gpu:')}
    if sum(free_types.values())<1:continue
    available.append(dict(node=d['NodeName'],partitions=d['Partitions'],state=d['State'],
        free_cpus=int(d['CPUTot'])-int(d['CPUAlloc']),
        free_mem_gib=(int(d['RealMemory'])-int(d['AllocMem']))/1024,
        os_free_mem_gib=int(d['FreeMem'])/1024 if d.get('FreeMem','').isdigit() else None,
        free_gpus=free_types))
(root/'resources.json').write_text(json.dumps(available,indent=2)+'\n')
print('AUDIT_ROOT='+str(root))
print((root/'baseline.tsv').read_text())
print('SOURCE_DIFFS='+json.dumps(diffs))
print('AVAILABLE='+json.dumps(available))
