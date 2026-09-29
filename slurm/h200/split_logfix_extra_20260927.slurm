#!/bin/bash
#SBATCH --job-name=ps-split-logfix
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_split_logfix_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_split_logfix_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 - <<'PY'
import csv, hashlib, json, os, subprocess
from pathlib import Path
assert os.environ['SLURM_JOB_ID']
base = Path('/usr/xtmp/lz280/packstar_pro6000_logfix_20260927')
original = base / 'postcheck_12706086'
package = base / 'incremental_package_12706086'
build = base / 'complete_12706035'
state = subprocess.check_output(['squeue','-h','-j','12706087','-o','%T'], text=True).strip()
assert state == 'PENDING', state
rows = list(csv.DictReader((original / 'new_cases.tsv').open(), delimiter='\t'))
assert [r['index'] for r in rows] == ['132', '134']
children = []
for row in rows:
    root = base / ('postcheck_split_' + os.environ['SLURM_JOB_ID'] + '_' + row['index'])
    root.mkdir()
    with (root / 'new_cases.tsv').open('w') as out:
        writer = csv.DictWriter(out, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    digest = hashlib.sha256((root / 'new_cases.tsv').read_bytes()).hexdigest()
    (root / 'cases.sha256').write_text(digest + '  new_cases.tsv\n')
    command = ['sbatch', '--parsable', '--hold', '--account=grisman',
               '--job-name=ps-logfix-' + row['index'],
               '--export=ALL,BUILD_ROOT=' + str(build) + ',FOLLOWUP_ROOT=' + str(root),
               str(package / 'node_log10_followup_20260927.slurm')]
    (root / 'submission_command.json').write_text(json.dumps(command, indent=2) + '\n')
    job = subprocess.check_output(command, text=True).strip().split(';')[0]
    assert job.isdigit()
    (root / 'submitted_job.txt').write_text(job + '\n')
    children.append(job)
    print(row['index'], row['design_id'], job, flush=True)
# Both replacements are registered before retiring the original queue entry.
state = subprocess.check_output(['squeue','-h','-j','12706087','-o','%T'], text=True).strip()
assert state == 'PENDING', state
subprocess.run(['scancel', '12706087'], check=True)
(original / 'submitted_job.txt').rename(original / 'superseded_job.txt')
(original / 'replacement_jobs.json').write_text(json.dumps(children) + '\n')
for job in children:
    subprocess.run(['scontrol', 'release', job], check=True)
print('SPLIT_COMPLETE', *children, flush=True)
PY
