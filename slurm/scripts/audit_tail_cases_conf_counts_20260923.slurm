#!/bin/bash
#SBATCH --job-name=tail_conf_counts
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
import json
import math
import os
import re
from pathlib import Path

root = Path('/usr/xtmp/lz280/packstar_pair38_20260921/launch_12677939/runs')
cases = [
    ('3cal', 'TRP', 'Complex', 'a0ead732c412f73f'),
    ('4wyu', 'ASN MET', 'Complex', '1a0fe07ba3d6f7db'),
    ('1gwc', 'ARG', 'Complex', 'b98b5d31c4ec0037'),
    ('3cal', 'TRP', 'Protein', 'e7ef13207bdd8a21'),
    ('5d68', 'ILE', 'Complex', 'e41f96d22045167'),
]
result = []
for design, sequence, role, identity in cases:
    run, = root.glob(design+'_full_pair-only_s42_*')
    lines = (run/'run.log').read_text().splitlines()
    stop, = [i for i,line in enumerate(lines) if 'identityHash='+identity in line]
    start = max(i for i in range(stop) if 'DP tables computed in' in lines[i])
    block = lines[start:stop+1]
    assert any('estimator activated for state='+role+',' in x for x in block)
    num_pos, = [int(re.search(r'isolated positions=\d+/(\d+)', x).group(1))
                for x in block if 'isolated positions=' in x]
    domains = {}
    for line in block:
        if 'Edge: lambda=' not in line and 'isolated positions=' not in line:
            continue
        for pos, count in re.findall(r'(\d+)\((\d+)rc\)', line):
            pos, count = int(pos), int(count)
            assert pos not in domains or domains[pos] == count
            domains[pos] = count
    assert set(domains) == set(range(num_pos)), (design, domains, num_pos)
    counts = [domains[p] for p in range(num_pos)]
    manifest = json.loads((run/'manifest.json').read_text())
    assert manifest['properties']['branchdp.cutoff.residualBudget'] == '1'
    pruning = [x for x in lines if 'Sparse emat stats:' in x]
    assert pruning and all('prunedRcPairs=0' in x for x in pruning)
    samples_path = run/'adaptive_frequency_severity'/('state-'+identity)/'frequency_severity_final_samples.tsv'
    with samples_path.open() as stream:
        samples = list(csv.DictReader(stream, delimiter='\t'))
    assert all(len(x['conf'].split(','))==num_pos for x in samples)
    row = dict(design=design, sequence=sequence, role=role, positions=num_pos,
               rc_counts=counts, conformations=math.prod(counts),
               final_draws=len(samples), final_unique_conformations=len({x['conf'] for x in samples}),
               final_tail_draws=sum(x['tailEvent']=='true' for x in samples),
               final_tail_unique_conformations=len({x['conf'] for x in samples if x['tailEvent']=='true'}),
               log=str(run/'run.log'), first_line=start+1, last_line=stop+1,
               state='state-'+identity)
    result.append(row)
    print(json.dumps(row), flush=True)
out = Path('/usr/xtmp/lz280') / ('packstar_tail_conf_counts_20260923_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
(out/'counts.json').write_text(json.dumps(result,indent=2)+'\n')
with (out/'counts.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(result[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(result)
report = ['# Conformation counts for five tail-exceeding final states', '',
          'Total = product of the per-position RC counts used in the sequence-specific DP. This counts discrete assignments; it is not a count of positive-weight or independently sampled conformations. All run logs report zero sparse cut-pair pruning.', '',
          '| System | Sequence | State | Positions | Total conformations | Final draws | Unique final conformations | Tail draws | Unique tail conformations |',
          '|---|---|---|---:|---:|---:|---:|---:|---:|']
for r in result:
    report.append('| '+' | '.join(str(r[k]) for k in ('design','sequence','role','positions','conformations','final_draws','final_unique_conformations','final_tail_draws','final_tail_unique_conformations'))+' |')
    report += []
report += ['', 'RC counts and exact source log line numbers are saved in counts.json and counts.tsv.']
(out/'report.md').write_text('\n'.join(report)+'\n')
print('OUTPUT='+str(out), flush=True)
PY
