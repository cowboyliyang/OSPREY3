"""Slurm-only historical workload and current scheduler inventory."""
import csv
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get('SLURM_JOB_ID')
root = Path('/usr/xtmp/lz280/packstar_pair38_rb24_20260922') / ('audit_' + os.environ['SLURM_JOB_ID'])
root.mkdir(parents=True, exist_ok=False)
history = []
for rb in (1, 2, 4):
    baseline = Path('/usr/xtmp/lz280/packstar_rb38_20260918/launch_12637884') / f'rb{rb}'
    for path in sorted(baseline.glob('*/manifest.json')):
        m = json.loads(path.read_text())
        design = m['design']['design_id']
        wall = (path.parent / 'wall.time').read_text()
        match = re.search(r'elapsed=([\d.]+) maxRssKiB=(\d+)', wall)
        assert m['options']['seed'] == 42, path
        assert m['properties']['branchdp.cutoff.residualBudget'] == str(rb)
        assert m['properties']['packstar.pac.frequencySeverity.tripleEta'] == 'true'
        output = path.parent / (design + '_packstar.csv')
        rr = list(csv.DictReader(output.open())) if output.exists() else []
        history.append(dict(rb=rb, design=design, index=m['options']['index'],
            elapsed_s=float(match[1]) if match else m.get('end_epoch', m['start_epoch'])-m['start_epoch'],
            peak_rss_gib=int(match[2])/2**20 if match else 0,
            wall_complete=bool(match),
            status=m['status'], exit_code=m.get('exit_code'), rows=len(rr),
            estimated=sum(all(r[s+'_status']=='Estimated' for s in ('prot','lig','comp')) for r in rr),
            path=str(path.parent)))
assert len(history) == 114
(root / 'history.json').write_text(json.dumps(history, indent=2) + '\n')
nodes = subprocess.check_output(['scontrol', 'show', 'nodes', '--oneliner'], text=True)
(root / 'nodes.txt').write_text(nodes)
available = []
for line in nodes.splitlines():
    d = dict(re.findall(r'(\w+)=([^ ]*)', line))
    if not set(d.get('Partitions','').split(',')).intersection({'grisman','compsci-gpu'}):
        continue
    if any(s in d.get('State','') for s in ['DRAIN','DOWN','NOT_RESPONDING']):
        continue
    def tres(s):
        return dict(part.split('=',1) for part in s.split(',') if '=' in part)
    cfg, alloc = tres(d.get('CfgTRES','')), tres(d.get('AllocTRES',''))
    free = {key[9:]:int(value)-int(alloc.get(key,0)) for key,value in cfg.items() if key.startswith('gres/gpu:')}
    if not sum(free.values()):
        continue
    available.append(dict(node=d['NodeName'], partitions=d['Partitions'], state=d['State'],
        free_cpus=int(d['CPUTot'])-int(d['CPUAlloc']),
        free_mem_gib=(int(d['RealMemory'])-int(d['AllocMem']))/1024,
        os_free_mem_gib=int(d['FreeMem'])/1024 if d.get('FreeMem','').isdigit() else None,
        free_gpus=free))
(root / 'resources.json').write_text(json.dumps(available, indent=2)+'\n')
print('AUDIT_ROOT', root)
print('design rb1_seconds rb2_seconds rb4_seconds max_rss_gib rb2_rows rb4_rows')
for design in sorted({r['design'] for r in history}, key=lambda design:max(r['elapsed_s'] for r in history if r['design']==design), reverse=True):
    h = {r['rb']:r for r in history if r['design']==design}
    print(design, *[round(h[rb]['elapsed_s']) for rb in (1,2,4)], round(max(r['peak_rss_gib'] for r in h.values()),1), h[2]['rows'], h[4]['rows'])
print('RESOURCES')
for r in available:
    print(r['node'], r['state'], r['free_cpus'], round(r['free_mem_gib']), r['free_gpus'])
print(subprocess.check_output(['df','-h','/usr/xtmp/lz280'], text=True))
print(subprocess.check_output(['df','-i','/usr/xtmp/lz280'], text=True))
