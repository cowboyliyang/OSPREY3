#!/bin/bash
#SBATCH --job-name=tail_mark_overlap
#SBATCH --account=grisman
#SBATCH --partition=compsci,grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/usr/xtmp/lz280
python3 - <<'PY'
import csv
import hashlib
import json
import math
import os
from pathlib import Path

base = Path('/usr/xtmp/lz280')
root = base/'packstar_pair38_20260921/launch_12677939/runs'
old = base/'bench_comparison'
frozen = json.loads((base/'recomb_pair_rb_20260922/collect_12687019/input_manifest.json').read_text())
cases = [
    ('3cal','TRP','comp','a0ead732c412f73f'),
    ('4wyu','ASN MET','comp','1a0fe07ba3d6f7db'),
    ('1gwc','ARG','comp','b98b5d31c4ec0037'),
    ('3cal','TRP','prot','e7ef13207bdd8a21'),
    ('5d68','ILE','comp','e41f96d22045167'),
]
def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))
def spec(path):
    with path.open() as stream:
        return {r[0]:r for r in csv.reader(x for x in stream if x.strip() and not x.startswith('#'))}
current = spec(Path('/home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/baseline38.csv'))
historical = spec(old/'design_specs_prepped.csv')
inputs = {}
results = []
for design, sequence, state, identity in cases:
    assert current[design][1:2]+current[design][5:7] == historical[design][1:2]+historical[design][5:7]
    run, = root.glob(design+'_full_pair-only_s42_*')
    pp = run/(design+'_packstar.csv')
    mp = old/'results'/(design+'_markstar.csv')
    for path in (pp, mp):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == frozen[str(path)]['sha256'], path
        inputs[str(path)] = dict(sha256=digest, bytes=path.stat().st_size)
    p, = [r for r in records(pp) if r['sequence'].strip()==sequence]
    m, = [r for r in records(mp) if r['sequence'].strip()==sequence]
    assert p['method']=='packstar' and m['method']=='markstar'
    assert p[state+'_status']==m[state+'_status']=='Estimated'
    pl, pu, ml, mu = [float(r[state+'_qstar_'+bound+'_log10']) for r in (p,m) for bound in ('lb','ub')]
    assert all(math.isfinite(x) for x in (pl,pu,ml,mu)) and pl<=pu and ml<=mu
    interval = run/'adaptive_frequency_severity'/('state-'+identity)/'frequency_severity_final_interval.tsv'
    with interval.open() as stream:
        final = {r['key']:r['value'] for r in csv.DictReader(stream,delimiter='\t')}
    assert abs(float(final['logZLower'])/math.log(10)-pl)<1e-6
    assert abs(float(final['logZUpper'])/math.log(10)-pu)<1e-6
    overlap = max(pl,ml)<=min(pu,mu)
    direction = 'overlap' if overlap else ('PACK upper < MARK lower' if pu<ml else 'MARK upper < PACK lower')
    gap = 0 if overlap else max(ml-pu,pl-mu)
    results.append(dict(design=design, sequence=sequence, state=state,
        pack_lower_log10=pl,pack_upper_log10=pu,mark_lower_log10=ml,mark_upper_log10=mu,
        overlap=overlap,direction=direction,gap_log10=gap,nearest_bound_ratio=10**gap,
        pack_epsilon=float(p[state+'_eps']),mark_epsilon=float(m[state+'_eps']),
        pack_source=str(pp),mark_source=str(mp)))
out = base/('packstar_tail_markstar_comparison_20260924_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
with (out/'comparison.tsv').open('w') as stream:
    writer = csv.DictWriter(stream,fieldnames=list(results[0]),delimiter='\t')
    writer.writeheader(); writer.writerows(results)
summary = dict(cases=results, disjoint=sum(not r['overlap'] for r in results),
    definition='Compare state-specific log10 partition-function intervals, not Kstar score intervals.',
    provenance='Historical MARK CSVs and current PACK CSVs match the paper input hashes. Mutable/flexible specs match. Historical PDB content checksum unavailable; sparse versus full target differences are not removed.',
    inputs=inputs,job=os.environ['SLURM_JOB_ID'])
(out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
lines = ['# Five tail-exceeding states: PACK versus historical MARK intervals', '',
         summary['definition'], '', summary['provenance'], '',
         '| System | Sequence | State | PACK log10 Z interval | MARK log10 Z interval | Overlap | Nearest bound ratio |',
         '|---|---|---|---|---|---|---:|']
for r in results:
    lines.append(f"| {r['design']} | {r['sequence']} | {r['state']} | [{r['pack_lower_log10']:.6f}, {r['pack_upper_log10']:.6f}] | [{r['mark_lower_log10']:.6f}, {r['mark_upper_log10']:.6f}] | {r['overlap']} | {r['nearest_bound_ratio']:.6g} |")
(out/'report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,indent=2))
print('OUTPUT='+str(out))
PY
