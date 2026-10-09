"""Maintain pool-wide admission for explicitly registered pending MARK* jobs."""
import json
import os
from pathlib import Path
import re
import subprocess
import time
import fcntl

assert os.environ.get('SLURM_JOB_ID'), 'Run through Slurm.'
ROOT = Path('/usr/xtmp/lz280/markstar_free_backfill_20261008')
ROOT.mkdir(exist_ok=True)
lock = (ROOT / 'controller.lock').open('w')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
OUT = ROOT / os.environ['SLURM_JOB_ID']
OUT.mkdir(exist_ok=True)
configuration = Path(os.environ.get('MARKSTAR_BACKFILL_CONFIG',
    '/home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/free_backfill_jobs.json'))
config = json.loads(configuration.read_text())
LIMITS = config['limits']
assert LIMITS == {'compsci': 10, 'grisman': 5}
JOBS = {r['job']: r['pool'] for r in config['jobs']}
ADOPT_HOLD = {r['job'] for r in config['jobs'] if r.get('adopt_hold')}
assert len(JOBS) == len(config['jobs'])
assert set(JOBS.values()) <= set(LIMITS)
(OUT / 'configuration.json').write_text(json.dumps(config, indent=2) + '\n')

def run(*args):
    return subprocess.check_output(args, text=True, timeout=45).strip()

def info(job):
    return dict(re.findall(r'(?:^|\s)(\w+)=(\S+)', run('scontrol', 'show', 'job', job, '-o')))

def queue():
    rows = {}
    for line in run('squeue', '-h', '-r', '-u', 'lz280', '-o', '%i|%j|%T|%P|%C|%a|%r').splitlines():
        j, name, state, pool, cpus, account, reason = line.split('|')
        if name.startswith('markstar') and cpus == '64' and pool in LIMITS:
            assert account == 'grisman'
            rows[j] = dict(state=state, pool=pool, reason=reason)
    return rows

def event(action, **fields):
    record = dict(time=time.time(), action=action, **fields)
    with (OUT / 'events.jsonl').open('a') as stream:
        stream.write(json.dumps(record) + '\n')
    print(json.dumps(record), flush=True)

# Completed tasks can disappear from scontrol; only inspect live registered jobs.
live = queue()
before = {j: info(j) for j in JOBS if j in live}
(OUT / 'before.json').write_text(json.dumps(before, indent=2) + '\n')
for j, s in before.items():
    assert s['Account'] == 'grisman' and s['NumCPUs'] == '64'
    assert s['Partition'] == JOBS[j]
    assert s['JobState'] in ('PENDING', 'RUNNING', 'COMPLETED')
    if s['JobState'] == 'PENDING':
        assert s['Reason'] != 'JobHeldUser' or j in ADOPT_HOLD, 'Do not take over an unregistered user hold.'

waiting = set()
for j, previous in before.items():
    current = info(j)
    if current['JobState'] != 'PENDING':
        continue
    if current['Reason'] != 'JobHeldUser':
        run('scontrol', 'hold', j)
    assert info(j)['JobState'] == 'PENDING'
    waiting.add(j)
    run('scontrol', 'update', 'JobId=' + j, 'Dependency=')
    after = info(j)
    assert after['Dependency'] in ('(null)', '0')
    for key in ('Command', 'NumCPUs', 'MinMemoryNode', 'TimeLimit', 'Partition', 'Account', 'ExcNodeList'):
        assert after[key] == before[j][key], (j, key)
    event('pooled', job=j, pool=JOBS[j])

while waiting:
    try:
        rows = queue()
        # Count running, completing, and already admitted pending tasks as occupied.
        occupied = {p: sum(r['pool'] == p and j not in waiting for j, r in rows.items()) for p in LIMITS}
        for j in JOBS:
            if j not in waiting:
                continue
            if j not in rows or rows[j]['state'] != 'PENDING':
                waiting.remove(j)
                continue
            p = JOBS[j]
            if occupied[p] < LIMITS[p]:
                assert info(j)['Reason'] == 'JobHeldUser'
                run('scontrol', 'release', j)
                waiting.remove(j)
                occupied[p] += 1
                event('released', job=j, pool=p, occupied=occupied[p])
        status = dict(controller=os.environ['SLURM_JOB_ID'], updated=time.time(),
                      waiting=sorted(waiting), occupied=occupied, limits=LIMITS,
                      audit=str(OUT))
        (ROOT / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
    except (subprocess.SubprocessError, OSError) as exc:
        event('retry', error=str(exc))
    if waiting:
        time.sleep(config.get('poll_seconds', 30))
event('all_managed_jobs_admitted')
