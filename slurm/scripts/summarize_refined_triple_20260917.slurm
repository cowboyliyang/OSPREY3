#!/bin/bash
#SBATCH --job-name=ps_refined_triple_summary
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${REFINED_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A${REFINED_JOB}
export PYTHONPYCACHEPREFIX="$ROOT/summary_pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/summarize_refined.slurm"
python3 - "$REFINED_JOB" <<'PY'
import csv, json, re, sys
from pathlib import Path
base = Path('/usr/xtmp/lz280/packstar_frontier_fullnode_20260917')
jobs = ['12626683', sys.argv[1]]
report = []
tables = []
for job in jobs:
    run = base/f'A{job}'/'runs'/f'2xgy_A{job}_T0'
    with (run/'2xgy_packstar.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 39
    good = [r['sequence'] for r in rows if all(r[f'{s}_status'] == 'Estimated' for s in ('prot', 'lig', 'comp'))]
    wall = float(re.search(r'elapsed=([\d.]+)', (run/'wall.time').read_text())[1])
    report.append(dict(job=job, run=str(run), wall_s=wall, minutes=wall/60, all_estimated=len(good),
                       estimated_sequences=sorted(good)))
    tables.append({r['sequence']: {k: v for k, v in r.items() if k != 'total_time_s'} for r in rows})
result = dict(runs=report, speedup=report[0]['wall_s']/report[1]['wall_s'],
              exact_csv_except_time=tables[0] == tables[1],
              regressions=sorted(set(report[0]['estimated_sequences']) - set(report[1]['estimated_sequences'])),
              scope='Full CSV comparison only; per-state numeric validation is in the four-case gate.')
(base/f'A{jobs[1]}'/'refinement_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
assert result['exact_csv_except_time'], 'Full CSV changed; inspect before accepting this full-run result.'
PY
