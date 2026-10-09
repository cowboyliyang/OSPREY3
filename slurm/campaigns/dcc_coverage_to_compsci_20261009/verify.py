"""Verify scheduler receipts, refreshed 89-design plan and DCC launch exclusions."""
import csv
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Verify through Slurm.'
here = Path(__file__).resolve().parent
inv = here.parent / 'current_inventory_20261006'
focus = here.parent / 'deadline_focus_20261006'
record = json.loads((here / 'registry.json').read_text())
with (inv / 'current_plan.tsv').open() as stream:
    plan = {r['design']: r for r in csv.DictReader(stream, delimiter='\t')}
assert len(plan) == 89 and len({r['system'] for r in plan.values()}) == 38
checks = []
for entry in record['jobs']:
    row = plan[entry['design']]
    assert row['job'] == entry['job'] and row['destination'] == 'compsci'
    assert row['state'] in ('备用暂停', '排队', '运行', '完成'), row
    text = subprocess.check_output(['scontrol', 'show', 'job', entry['job'], '-o'], text=True)
    state = dict(re.findall(r'(?:^|\s)(\w+)=(\S+)', text))
    assert state['Account'] == 'grisman' and state['NumCPUs'] == '64' and state['MinMemoryNode'] == '192G'
    assert state['Partition'] == 'compsci' and state['TimeLimit'] == '14-00:00:00'
    assert state['Dependency'] == '(null)'
    for filename in ('dcc_new_cases.txt', 'dcc_keep_existing_pending_cases.txt'):
        assert entry['design'] not in (focus / filename).read_text().splitlines()
    for filename in ('dcc_launch_priority.tsv', 'dcc_new_designs.tsv'):
        with (focus / filename).open() as stream:
            assert entry['design'] not in {r['design_id'] for r in csv.DictReader(stream, delimiter='\t')}
    checks.append(dict(design=entry['design'], job=entry['job'], state=state['JobState'], reason=state['Reason']))
controller = subprocess.check_output(['scontrol', 'show', 'job', record['controller_job'], '-o'], text=True)
assert 'JobState=RUNNING' in controller or 'JobState=PENDING' in controller
summary = dict(status='PASS', plan_designs=89, systems=38, jobs=checks, controller_job=record['controller_job'],
               remote_scheduler_accessed=False, remote_state='UNKNOWN', preflights=record['verification'])
(here / 'verification.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
(Path(record['prep']) / 'verification.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
print('VERIFIED', json.dumps(summary, ensure_ascii=False), flush=True)
