"""Rewire five never-started MARK* tasks without changing the lane count."""
import csv
import json
import os
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

assert os.environ.get('SLURM_JOB_ID'), 'Run this management task through Slurm.'
campaign = Path(__file__).resolve().parent
repo_campaigns = campaign.parent
plan = json.loads((campaign / 'plan.json').read_text())
audit = Path('/usr/xtmp/lz280/markstar_coverage_lane_rebalance_20261006') / os.environ['SLURM_JOB_ID']
audit.mkdir(parents=True, exist_ok=False)
commands = []

def run(*args):
    result = subprocess.run(args, capture_output=True, text=True)
    commands.append({'arguments': list(args), 'exit_code': result.returncode,
                     'stdout': result.stdout, 'stderr': result.stderr})
    (audit / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    result.check_returncode()
    return result.stdout

def state(job):
    raw = run('scontrol', 'show', 'job', job, '--oneliner')
    return dict(re.findall(r'(?:^|\s)([A-Za-z][A-Za-z0-9_]*)=(\S+)', raw))

def parent(fields):
    found = re.fullmatch(r'afterany:([0-9]+_[0-9]+)(?:\(unfulfilled\))?', fields['Dependency'])
    assert found, fields['Dependency']
    return found[1]

def never_started(fields):
    assert fields['JobState'] == 'PENDING', fields
    assert fields['RunTime'] == '00:00:00' and fields.get('Restarts', '0') == '0', fields
    assert fields['StartTime'] in ('Unknown', 'N/A'), fields
    assert fields['Account'] == 'grisman' and fields['NumCPUs'] == '64', fields
    assert fields['Partition'] == 'compsci', fields

edge_path = repo_campaigns / 'deadline_focus_20261006/local_lane_dependencies.tsv'
with edge_path.open(newline='') as handle:
    edges = list(csv.DictReader(handle, delimiter='\t'))
edge_map = {r['job']: r['parent'] for r in edges}
for change in plan['changes']:
    assert edge_map[change['job']] == change['old_parent'], change
    edge_map[change['job']] = change['new_parent']
assert max(Counter(edge_map.values()).values()) == 1, 'A lane would branch.'
for child in edge_map:
    seen = set()
    node = child
    while node in edge_map:
        assert node not in seen, 'Dependency cycle.'
        seen.add(node)
        node = edge_map[node]

before = {c['job']: state(c['job']) for c in plan['changes']}
for change in plan['changes']:
    fields = before[change['job']]
    never_started(fields)
    assert parent(fields) == change['old_parent'], change
    assert fields['Reason'] != 'JobHeldUser', 'Do not release a pre-existing user hold.'
(audit / 'before.json').write_text(json.dumps(before, indent=2) + '\n')
held, updated = [], []
try:
    for change in plan['changes']:
        never_started(state(change['job']))
        run('scontrol', 'hold', change['job'])
        held.append(change['job'])
        never_started(state(change['job']))
    for change in plan['changes']:
        never_started(state(change['job']))
        run('scontrol', 'update', 'JobId=' + change['job'], 'Dependency=afterany:' + change['new_parent'])
        updated.append(change)
        assert parent(state(change['job'])) == change['new_parent']
except Exception:
    rollback_ok = True
    for change in reversed(updated):
        try:
            never_started(state(change['job']))
            run('scontrol', 'update', 'JobId=' + change['job'], 'Dependency=afterany:' + change['old_parent'])
        except Exception:
            rollback_ok = False
    if rollback_ok:
        for job in reversed(held):
            run('scontrol', 'release', job)
    (audit / 'failure.json').write_text(json.dumps({'rollback_succeeded': rollback_ok, 'held': held}, indent=2))
    raise

for job in reversed(held):
    run('scontrol', 'release', job)
after = {c['job']: state(c['job']) for c in plan['changes']}
for change in plan['changes']:
    assert parent(after[change['job']]) == change['new_parent']
    assert after[change['job']]['Reason'] != 'JobHeldUser'
(audit / 'after.json').write_text(json.dumps(after, indent=2) + '\n')

applied = datetime.now(timezone.utc).isoformat()
for row in edges:
    if row['job'] in {c['job'] for c in plan['changes']}:
        row['parent'] = edge_map[row['job']]
        row['action'] = 'COVERAGE_LANE_REBALANCE_20261006'
with edge_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['job', 'parent', 'action'], delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(edges)

change_map = {c['job']: c for c in plan['changes']}
submission_path = repo_campaigns / 'deadline_focus_20261006/cpu_continuation_submissions.tsv'
with submission_path.open(newline='') as handle:
    submissions = list(csv.DictReader(handle, delimiter='\t'))
for row in submissions:
    if row['job'] in change_map:
        row['after_job'] = change_map[row['job']]['new_parent']
    elif row['job'] == '12827032_2':
        row['after_job'] = '12827711_0'
with submission_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(submissions[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(submissions)

record = dict(plan, status='APPLIED', management_job=os.environ['SLURM_JOB_ID'],
              applied_utc=applied, audit=str(audit))
(audit / 'summary.json').write_text(json.dumps(record, indent=2) + '\n')
(campaign / 'registry.json').write_text(json.dumps(record, indent=2) + '\n')
for relative in ['deadline_focus_20261006/registry.json', 'deadline_focus_20261006/policy.json']:
    path = repo_campaigns / relative
    content = json.loads(path.read_text())
    lists = [content.get('new_local_jobs', []), content.get('cpu_continuation', {}).get('local_submitted', [])]
    for entries in lists:
        for entry in entries:
            if entry['job'] in change_map:
                entry.setdefault('original_after_job', entry['after_job'])
                entry['after_job'] = change_map[entry['job']]['new_parent']
                entry['dependency_updated_utc'] = applied
            elif entry['job'] == '12827032_2':
                entry.setdefault('original_after_job', entry['after_job'])
                entry['after_job'] = '12827711_0'
    content['coverage_lane_rebalance_registry'] = str(campaign / 'registry.json')
    if path.name == 'policy.json':
        content['two_tier_coverage_goal']['observed_tiers_due'] = '2026-10-16'
        content['two_tier_coverage_goal']['final_round_start_date'] = '2026-10-16'
    path.write_text(json.dumps(content, indent=2) + '\n')
path = repo_campaigns / 'two_tier_coverage_20261006/registry.json'
content = json.loads(path.read_text())
content.setdefault('original_after_job', content['after_job'])
content['after_job'] = '12812439_7'
content['dependency_updated_utc'] = applied
content['coverage_lane_rebalance_registry'] = str(campaign / 'registry.json')
path.write_text(json.dumps(content, indent=2) + '\n')
print(json.dumps(record, indent=2), flush=True)
