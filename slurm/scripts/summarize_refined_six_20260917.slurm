#!/bin/bash
#SBATCH --job-name=ps_refined_six_summary
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
: "${OTHER_JOB:?}"
: "${TWOXGY_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A${OTHER_JOB}
mkdir -p "$ROOT"
export PYTHONPYCACHEPREFIX="$ROOT/summary_pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/summarize_six.slurm"
python3 - "$OTHER_JOB" "$TWOXGY_JOB" <<'PY'
import csv, json, re, sys
from pathlib import Path
base = Path('/usr/xtmp/lz280/packstar_frontier_fullnode_20260917')
other, twoxgy = sys.argv[1:]
systems = ['2xgy', '4u3s', '1a0r', '3ma2', '4wyu', '5a6y']
results = []
for task, system in enumerate(systems):
    entries, tables = [], []
    for job in ['12626683', twoxgy if task == 0 else other]:
        run = base/f'A{job}'/'runs'/f'{system}_A{job}_T{task}'
        path = run/f'{system}_packstar.csv'
        rows = list(csv.DictReader(path.open())) if path.exists() else []
        good = sorted(r['sequence'] for r in rows if all(r[f'{s}_status'] == 'Estimated' for s in ('prot', 'lig', 'comp')))
        wallfile = run/'wall.time'
        match = re.search(r'elapsed=([\d.]+)', wallfile.read_text()) if wallfile.exists() else None
        manifest = run/'run_manifest.tsv'
        metadata = dict(line.split('\t', 1) for line in manifest.read_text().splitlines()) if manifest.exists() else {}
        wall = float(match[1]) if match else None
        entries.append(dict(job=job, run=str(run), wall_s=wall,
                            minutes=wall/60 if wall else None, rows=len(rows),
                            all_estimated=len(good), estimated_sequences=good,
                            java_status=metadata.get('javaStatus', 'missing')))
        tables.append({r['sequence']: {k: v for k, v in r.items() if k != 'total_time_s'} for r in rows})
    old, new = entries
    differences = []
    for seq in sorted(set(tables[0]) | set(tables[1])):
        a, b = tables[0].get(seq, {}), tables[1].get(seq, {})
        fields = {k: [a.get(k), b.get(k)] for k in set(a) | set(b) if a.get(k) != b.get(k)}
        if fields: differences.append(dict(sequence=seq, fields=fields))
    results.append(dict(system=system, baseline=old, optimized=new,
                        wall_ratio=old['wall_s']/new['wall_s'] if old['wall_s'] and new['wall_s'] else None,
                        exact_csv_except_time=bool(tables[1]) and not differences,
                        lost_successes=sorted(set(old['estimated_sequences'])-set(new['estimated_sequences'])),
                        gained_successes=sorted(set(new['estimated_sequences'])-set(old['estimated_sequences'])),
                        csv_differences=differences))
root = base/f'A{other}'
(root/'six_system_comparison.json').write_text(json.dumps(results, indent=2)+'\n')
lines = ['Latest optimized full-node runs versus build 12626676.', '',
         '| System | Previous min | Latest min | Latest success | Wall ratio | CSV equal except time |',
         '|---|---:|---:|---:|---:|---|']
def fmt(value): return f'{value:.2f}' if value is not None else 'missing'
for r in results:
    a, b = r['baseline'], r['optimized']
    lines.append(f"| {r['system']} | {fmt(a['minutes'])} | {fmt(b['minutes'])} | {b['all_estimated']}/{b['rows']} | {fmt(r['wall_ratio'])} | {r['exact_csv_except_time']} |")
lines += ['', 'Wall ratios compare whole runs; inspect success sets and CSV differences before interpreting them.',
          'CSV equality is a formatted full-run check; per-state numerical validation is recorded separately in gate 12626763.']
(root/'six_system_comparison.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
PY
