"""Read-only input analysis; generate review drafts only, under Slurm."""
import csv, io, json, os, subprocess
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
repo = Path('/home/users/lz280/OSPREY3-fresh-packstar')
out = Path('/usr/xtmp/lz280') / ('frontier_review_' + os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
def read(path):
    return list(csv.DictReader(io.StringIO(path.read_text().replace('\\t', '\t')), delimiter='\t'))
screen = {r['design']:r for r in read(Path('/usr/xtmp/lz280/packstar_frontier_screen_20260917/12620972/all38.tsv'))}
existing = sum((read(repo/'slurm/h200'/f) for f in ['frontier.tsv','frontier_add4.tsv','frontier_small4.tsv']), [])
covered = {r['system'] for r in existing}
specs = {r[0]:r for r in csv.reader((repo/'slurm/h200/baseline38.csv').open()) if r and not r[0].startswith('#')}
remaining = sorted(set(specs)-covered)
assert len(remaining)==18
rows=[]
print('| 系统 | 原始点数 | 历史小时 | 低档 点/小时 | 中档 点/小时 | 高档 点/小时 | 高档delta |')
print('|---|---:|---:|---:|---:|---:|---:|')
for s in remaining:
    spec=specs[s]; mut=spec[5].split(';'); flex=list(dict.fromkeys(spec[6].split(';')))
    n=len(mut)+len(flex); h=float(screen[s]['mark_seconds'])/3600
    k=0
    while h*3**k <= 24: k+=1
    # If original already exceeds a day, shrink two protein-side WT sites.
    levels=[k-2,k-1,k]
    scans={int(r['k']):r for r in read(Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config/expansion_scans')/f'{s}_scan.tsv')}
    cells=[]
    for tier,d in zip(['low','middle','high'],levels):
        removal=[r for r in reversed(flex) if r[0] in {m[0] for m in mut}][:-d] if d<0 else []
        added=scans[d]['added_residues'].split(';') if d>0 else []
        chosen=[r for r in flex if r not in removal]+added
        assert len(chosen)==len(flex)+d and len(set(chosen+mut))==n+d
        row=dict(review_index=len(rows),design_id=f'{s}_flex_'+('m'+str(-d) if d<0 else 'p'+str(d)),system=s,tier=tier,flex_delta=d,total_positions=n+d,mutable=';'.join(mut),flexible=';'.join(chosen),hours_3x=h*3**d,hours_2x=h*2**d,partition='compsci,grisman' if tier=='high' else 'grisman',node_rule='exclude fennario-[01-06]' if tier=='high' else 'only fennario-[01-06]',status='DRAFT_NOT_PREFLIGHTED',job_id='')
        rows.append(row); cells.append(f'{n+d} / {h*3**d:.2f}')
    print(f'| {s} | {n} | {h:.2f} | '+ ' | '.join(cells)+f' | {k:+d} |')
for s,h in [('4znc',5.63),('3gxu',17.98),('4wem',18.30),('2rfd',23.42)]:
    old=max((r for r in existing if r['system']==s),key=lambda r:int(r['flex_delta']))
    d=int(old['flex_delta'])+1
    scans={int(r['k']):r for r in read(Path('/usr/xtmp/lz280/frontier_dcc_20260805_v2/config/expansion_scans')/f'{s}_scan.tsv')}
    prev=scans[d-1]['added_residues'].split(';'); added=scans[d]['added_residues'].split(';')
    extra=[r for r in added if r not in prev]; assert len(extra)==1
    flex=old['flexible'].split(';')+extra
    rows.append(dict(review_index=len(rows),design_id=f'{s}_flex_p{d}',system=s,tier='extra_plus1',flex_delta=d,total_positions=len(flex)+len(old['mutable'].split(';')),mutable=old['mutable'],flexible=';'.join(flex),hours_3x=h*3,hours_2x=h*2,partition='compsci',node_rule='compsci only',status='DRAFT_NOT_PREFLIGHTED',job_id=''))
    print('EXTRA',json.dumps(rows[-1]))
with (out/'draft_designs.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
nodes=subprocess.check_output(['scontrol','show','nodes','-o'],text=True)
(out/'nodes.txt').write_text(nodes)
totals={}
for line in nodes.splitlines():
    a=dict(x.split('=',1) for x in line.split() if '=' in x)
    if not set(a.get('Partitions','').split(','))&{'compsci','grisman'}:continue
    group='fennario' if a['NodeName'].startswith('fennario') else a['Partitions']
    slots=max(0,min((int(a['CPUEfctv'])-int(a['CPUAlloc']))//16,(int(a['RealMemory'])-int(a['AllocMem']))//98304))
    if any(x in a['State'] for x in ['DRAIN','DOWN','NOT_RESPONDING']):slots=0
    totals[group]=totals.get(group,0)+slots
    print('CAPACITY',a['NodeName'],slots,a['State'])
print('TOTAL_CAPACITY',totals)
print('OUTPUT',out)
