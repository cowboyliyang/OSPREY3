#!/bin/bash
#SBATCH --job-name=compare_4wyu_z
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
python3 - <<'PY'
import csv,json,math,statistics
from pathlib import Path
ref=Path('/usr/xtmp/lz280/packstar_million_census_20260918/run_12633933/merged_A12634210')
pack=Path('/usr/xtmp/lz280/packstar_million_coverage_20260918/run_12634011')
mark=Path('/usr/xtmp/lz280/markstar_million5_20260918/build_12638085/runs/4wyu_J12638087')
assert (ref/'READY').exists() and (mark/'READY').exists()
truth=json.loads((ref/'reference.json').read_text())
coverage=json.loads((pack/'coverage.json').read_text())
assert (mark/'layout.tsv').read_bytes()==(ref/'layout.tsv').read_bytes()
config=next(csv.DictReader((mark/'config.tsv').open(),delimiter='\t'))
assert config==truth['config']
m=next(csv.DictReader((mark/'result.tsv').open(),delimiter='\t'))
z=truth['logZ']; lo=float(m['logLower']); hi=float(m['logUpper'])
times=[]; widths=[]
for seed in range(10000,10100):
    p=pack/f'seed_{seed}'
    assert (p/'layout.tsv').read_bytes()==(ref/'layout.tsv').read_bytes()
    r=next(csv.DictReader((p/'result.tsv').open(),delimiter='\t'))
    widths.append(-math.expm1(float(r['logLower'])-float(r['logUpper'])))
    fields=dict(item.split('=',1) for item in (p/'wall.time').read_text().split() if '=' in item)
    times.append(float(fields['elapsed']))
mt=dict(item.split('=',1) for item in (mark/'wall.time').read_text().split() if '=' in item)
report=dict(system='4wyu',count=truth['total'],logZ=z,Z=math.exp(z),
    markstar=dict(logLower=lo,logUpper=hi,covered=lo<=z<=hi,
        epsilon=-math.expm1(lo-hi),process_seconds=float(mt['elapsed']),
        lower_over_Z=math.exp(lo-z),upper_over_Z=math.exp(hi-z)),
    packstar=dict(coverage=coverage['counts'],wilson95=coverage['sampled_coverage_wilson95'],
        process_seconds_min=min(times),process_seconds_median=statistics.median(times),
        process_seconds_max=max(times),epsilon_median=statistics.median(widths)),
    pack_median_over_mark_time=statistics.median(times)/float(mt['elapsed']))
out=pack/'three_way_comparison.json'
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
PY
