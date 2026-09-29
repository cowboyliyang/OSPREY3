"""After the original workers finish, submit only cases missed at primary launch."""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

assert os.environ['SLURM_JOB_ID']
root = Path(os.environ['FOLLOWUP_ROOT'])
build = Path('/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/complete_12706035')
primary = Path('/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/node_12706036')
package = Path(__file__).parent

def finish(status, **fields):
    result = dict(status=status, primary_job='12706036', **fields)
    (root / 'decision.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)

for script in ('audit_all_log10_20260927.py', 'select_all_log10_reruns_20260927.py'):
    args = [sys.executable, str(build / 'source/slurm/h200' / script), str(root)]
    if script.startswith('select_'): args.append(str(build))
    with (root / (script + '.log')).open('w') as stream:
        subprocess.run(args, stdout=stream, stderr=subprocess.STDOUT, check=True)
with (root / 'cases.tsv').open() as stream:
    reader = csv.DictReader(stream, delimiter='\t')
    fields = reader.fieldnames
    candidates = list(reader)

query = subprocess.run(['squeue', '-h', '-j', '12706036', '-o', '%T'], text=True, capture_output=True)
state = query.stdout.strip() if query.returncode == 0 else ''
if not state:
    # squeue can reject an ID after a completed job leaves its live records.
    # Accounting retains the terminal state needed for the deduplication audit.
    states = subprocess.check_output(['sacct', '-X', '-n', '-P', '-j', '12706036', '--format=State'], text=True).strip().splitlines()
    state = states[0].split('|')[0].split()[0] if states else 'UNKNOWN'

if state == 'PENDING':
    # The independently queued primary worker scans all final original results
    # at launch. Submitting another GPU worker here would duplicate that work.
    finish('PRIMARY_WILL_INCLUDE_FINAL_SCAN', candidate_count=len(candidates), primary_state=state)
    raise SystemExit(0)
if state not in ('RUNNING', 'COMPLETING', 'COMPLETED'):
    finish('PRIMARY_STATE_REQUIRES_REVIEW', candidate_count=len(candidates), primary_state=state)
    raise SystemExit(1)

for attempt in range(12):
    if (primary / 'cases.tsv').exists() and (primary / 'status.tsv').exists(): break
    time.sleep(5)
else:
    finish('PRIMARY_CASE_LIST_UNAVAILABLE', candidate_count=len(candidates), primary_state=state)
    raise SystemExit(1)

covered = {r['index'] for r in csv.DictReader((primary / 'cases.tsv').open(), delimiter='\t')}
# Include supplementary workers already submitted by an earlier audit.
for prior in sorted(root.parent.glob('postcheck_*/submitted_job.txt')):
    prior_cases = prior.parent / 'new_cases.tsv'
    job = prior.read_text().strip()
    states = subprocess.check_output(['sacct', '-X', '-n', '-P', '-j', job,
                                      '--format=State'], text=True).strip().splitlines()
    prior_state = states[0].split('|')[0].split()[0] if states else 'UNKNOWN'
    if prior_state not in ('PENDING', 'RUNNING', 'COMPLETING', 'COMPLETED'):
        raise RuntimeError(f'Supplementary job {job} needs review: {prior_state}')
    covered.update(r['index'] for r in csv.DictReader(prior_cases.open(), delimiter='\t'))

# Record completed-output cases added since the primary launch audit.
previous = {r['index']: r for r in csv.DictReader((primary / 'all_162.tsv').open(), delimiter='\t')}
current = list(csv.DictReader((root / 'all_162.tsv').open(), delimiter='\t'))
delta = [r for r in current if r['stage'] == 'CSV_WRITTEN'
         and previous[r['index']]['stage'] != 'CSV_WRITTEN']
with (root / 'newly_written.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(delta)
for r in delta:
    print('NEW_CSV', r['index'], r['design_id'], 'gap_s=' + r['output_gap_s'],
          'extreme=' + r['any_extreme'], flush=True)
new = [r for r in candidates if r['index'] not in covered]
with (root / 'new_cases.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(new)
if not new:
    finish('NO_ADDITIONAL_CASES', candidate_count=len(candidates), already_covered=len(covered))
    raise SystemExit(0)

checksum = hashlib.sha256((root / 'new_cases.tsv').read_bytes()).hexdigest()
(root / 'cases.sha256').write_text(checksum + '  new_cases.tsv\n')
command = ['sbatch', '--parsable', '--account=grisman',
           '--export=ALL,BUILD_ROOT=' + str(build) + ',FOLLOWUP_ROOT=' + str(root),
           str(package / 'node_log10_followup_20260927.slurm')]
(root / 'submission_command.json').write_text(json.dumps(command, indent=2) + '\n')
child = subprocess.check_output(command, text=True).strip().split(';')[0]
assert child.isdigit(), child
(root / 'submitted_job.txt').write_text(child + '\n')
finish('ADDITIONAL_CASES_SUBMITTED', candidate_count=len(candidates), already_covered=len(covered),
       additional_case_count=len(new), submitted_job=child)
