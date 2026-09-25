#!/bin/bash
#SBATCH --job-name=ps_gpu_ccd_summary
#SBATCH --account=users
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${GPU_CCD_JOB:?}"
ROOT=/usr/xtmp/lz280/packstar_frontier_gpu_ccd_20260917/A${GPU_CCD_JOB}
mkdir -p "$ROOT"
export PYTHONPYCACHEPREFIX="$ROOT/summary_pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/summarize_gpu_ccd.slurm"
python3 - "$GPU_CCD_JOB" <<'PY'
import csv, json, math, re, sys
from pathlib import Path
job = sys.argv[1]
base = Path('/usr/xtmp/lz280/packstar_frontier_fullnode_20260917')
root = Path('/usr/xtmp/lz280/packstar_frontier_gpu_ccd_20260917')/f'A{job}'
report = []
def read(run):
    paths = list(run.glob('*_packstar.csv'))
    rows = list(csv.DictReader(paths[0].open())) if paths else []
    tables = {r['sequence']: {k:v for k,v in r.items() if k != 'total_time_s'} for r in rows}
    good = sorted(r['sequence'] for r in rows if all(r[f'{s}_status'] == 'Estimated' for s in ('prot','lig','comp')))
    wall = run/'wall.time'
    match = re.search(r'elapsed=([\d.]+)', wall.read_text()) if wall.exists() else None
    mf = run/'run_manifest.tsv'
    meta = dict(line.split('\t',1) for line in mf.read_text().splitlines()) if mf.exists() else {}
    return dict(run=str(run), wall_s=float(match[1]) if match else None,
                rows=len(rows), successes=good, java_status=meta.get('javaStatus','missing')), tables
for task, system in enumerate(['2xgy','4u3s','1a0r','3ma2','4wyu','5a6y']):
    cpu_job = '12626765' if task == 0 else '12627258'
    cpu, a = read(base/f'A{cpu_job}'/'runs'/f'{system}_A{cpu_job}_T{task}')
    gpu, b = read(root/'runs'/f'{system}_A{job}_T{task}')
    diffs, numeric_max = [], {}
    for seq in sorted(set(a)|set(b)):
        fields = {}
        for k in sorted(set(a.get(seq,{}))|set(b.get(seq,{}))):
            old, new = a.get(seq,{}).get(k), b.get(seq,{}).get(k)
            if old == new: continue
            fields[k] = [old,new]
            try:
                x,y = float(old),float(new)
                if math.isfinite(x) and math.isfinite(y):
                    numeric_max[k] = max(numeric_max.get(k,0),abs(x-y))
            except (TypeError,ValueError): pass
        if fields: diffs.append(dict(sequence=seq,fields=fields))
    report.append(dict(system=system,cpu=cpu,gpu=gpu,
        complete=gpu['java_status']=='0' and gpu['rows']==cpu['rows'],
        speedup=cpu['wall_s']/gpu['wall_s'] if cpu['wall_s'] and gpu['wall_s'] else None,
        exact_csv_except_time=bool(b) and a==b,
        lost_successes=sorted(set(cpu['successes'])-set(gpu['successes'])),
        gained_successes=sorted(set(gpu['successes'])-set(cpu['successes'])),
        max_abs_numeric_differences=numeric_max,csv_differences=diffs))
(root/'cpu_gpu_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['FP64 GPU CCD versus latest FP64 CPU CCD, same build 12626762 and full-node allocation.','',
       '| System | CPU min | GPU min | CPU success | GPU success | Speedup | Complete |',
       '|---|---:|---:|---:|---:|---:|---|']
def minutes(x): return f'{x/60:.2f}' if x is not None else 'missing'
for r in report:
    a,b=r['cpu'],r['gpu']
    ratio=f"{r['speedup']:.3f}" if r['speedup'] is not None else 'missing'
    lines.append(f"| {r['system']} | {minutes(a['wall_s'])} | {minutes(b['wall_s'])} | {len(a['successes'])}/{a['rows']} | {len(b['successes'])}/{b['rows']} | {ratio} | {r['complete']} |")
lines += ['', 'Compare success sets and numerical differences in the JSON before interpreting speedups. GPU and CPU minimization can differ despite both using FP64.']
(root/'cpu_gpu_comparison.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
PY
