#!/bin/bash
#SBATCH --job-name=ps_rb1_baseline_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 - <<'PY'
import csv, json, os
from pathlib import Path
base = Path('/usr/xtmp/lz280/packstar_rb38_20260918/launch_12637884/rb1')
systems = []
for path in sorted(base.glob('*/manifest.json')):
    meta = json.loads(path.read_text())
    design = meta['design']['design_id']
    result = path.parent / (design + '_packstar.csv')
    rows = list(csv.DictReader(result.open())) if result.exists() else []
    expected = max(meta['expected_rows'], len(rows))
    success = sum(all(r[s+'_status']=='Estimated' for s in ('prot','lig','comp')) for r in rows)
    systems.append(dict(design=design, expected=expected, observed=len(rows), estimated=success,
                        all_estimated=len(rows)==expected==success, process_exit=meta.get('exit_code')))
summary = dict(systems=len(systems), all_estimated_systems=sum(r['all_estimated'] for r in systems),
    expected_sequences=sum(r['expected'] for r in systems), observed_sequences=sum(r['observed'] for r in systems),
    estimated_sequences=sum(r['estimated'] for r in systems),
    fully_successful_systems=[r['design'] for r in systems if r['all_estimated']],
    systems_detail=systems)
root=Path('/usr/xtmp/lz280/packstar_pair38_20260921')
(root / ('baseline_summary_'+os.environ['SLURM_JOB_ID']+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
PY
