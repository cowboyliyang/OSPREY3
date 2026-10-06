"""Audit the accepted handoff and plan (do not launch) the 24 MARK* cases."""
import collections
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import time

assert os.environ.get('SLURM_JOB_ID')
here = Path(__file__).parent
repo = here.parents[2]
out = Path('/usr/xtmp/lz280/frontier_extension7_20261006')/('markstar_plan_'+os.environ['SLURM_JOB_ID'])
out.mkdir()
campaign = repo/'slurm/campaigns/markstar_dcc_20261005'
handoff = json.loads((campaign/'rebalance_20261006.json').read_text())
queue_raw = subprocess.check_output(['squeue','-h','-r','-j','12812439','-o','%i|%T|%N|%r'], text=True)
(out/'queue_before.txt').write_text(queue_raw)
queue = {p[0]:p[1:] for p in (line.split('|') for line in queue_raw.splitlines())}
ids = [r['local_job'] for r in handoff['additions']]
assert set(ids)=={f'12812439_{i}' for i in [14,15,18,19,23,27,33,37,44]}
assert all(queue[j][0]=='PENDING' and queue[j][2]=='JobHeldUser' for j in ids)
running_before = {j:r for j,r in queue.items() if r[0]=='RUNNING'}
# Explicit user report on 2026-10-06 confirms external submission of all nine.
# Cancel only the previously held duplicate copies; every running job is excluded.
subprocess.run(['scancel', *ids], check=True)
after_raw = subprocess.check_output(['squeue','-h','-r','-j','12812439','-o','%i|%T|%N|%r'], text=True)
(out/'queue_after.txt').write_text(after_raw)
after = {p[0]:p[1:] for p in (line.split('|') for line in after_raw.splitlines())}
assert all(j not in after or after[j][0]=='CANCELLED' for j in ids)
assert all(j in after and after[j][0:2]==r[0:2] for j,r in running_before.items())
handoff['external_acceptance'] = dict(source='User confirmed submission on 2026-10-06; no new completed designs',
                                      external_job_ids='Not supplied', reported_completed_count=24,
                                      local_duplicate_cleanup_job=os.environ['SLURM_JOB_ID'], audit=str(out))
for r in handoff['additions']:
    r.update(local_status='CANCELLED_DUPLICATE', external_status='SUBMITTED_USER_CONFIRMED')
(campaign/'rebalance_20261006.json').write_text(json.dumps(handoff, indent=2)+'\n')
readme = (campaign/'README.md').read_text()
old = ('The nine local tasks are held for handoff so they do not start duplicate\n'
       'work. They have not been cancelled. After DCC acceptance is confirmed, retire\n'
       'the held local copies. Keep newly transferred tiers of each system on one CPU\n'
       'model and record actual CPU topology and binding.')
new = ('The user confirmed DCC submission of all nine additions on October 6, with\n'
       'no newly completed designs. The nine held local duplicate tasks were then\n'
       'cancelled after verifying that each was still pending and held. Running\n'
       'workloads were preserved. External job IDs have not been supplied. Keep\n'
       'newly transferred tiers of each system on one CPU model and record actual\n'
       'CPU topology and binding.')
assert readme.count(old)==1
(campaign/'README.md').write_text(readme.replace(old,new))

inventory = json.loads(Path('/usr/xtmp/lz280/markstar_extension_inventory_12826374/inventory.json').read_text())
policy = json.loads(Path('/home/users/lz280/notes/markstar_hardware_policy_12812439.json').read_text())
tasks = {r['design']:r for r in policy['tasks']}
new_designs = list(csv.DictReader((here/'new_designs.tsv').open(), delimiter='\t'))
assert len(new_designs)==7 and len(inventory['remaining'])==17
existing_local = set((campaign/'local_cases.txt').read_text().splitlines())
existing_dcc = set((campaign/'dcc_cases.txt').read_text().splitlines())
assert len(existing_local)==48 and len(existing_dcc)==91
active = {}
for path in Path('/usr/xtmp/lz280/markstar_local57_cpu64_20261005/A12812439').glob('*/run_manifest.json'):
    m=json.loads(path.read_text())
    active[m['design_id']]=dict(node=m['node'], status=m['status'],
                              group='fennario' if m['node'].startswith('fennario') else 'compsci')
nodes_raw = subprocess.check_output(['scontrol','show','nodes','-o'], text=True)
(out/'nodes.txt').write_text(nodes_raw)
slots = collections.Counter()
node_slots = []
for line in nodes_raw.splitlines():
    n=dict(t.split('=',1) for t in line.split() if '=' in t)
    name=n['NodeName']
    if not (re.fullmatch(r'fennario-0[1-6]',name) or re.fullmatch(r'compsci-cluster-fitz-(3[5-9]|4[0-4])',name)):
        continue
    if any(s in n['State'] for s in ['DOWN','DRAIN']):
        continue
    group='fennario' if name.startswith('fennario') else 'compsci'
    own=sum(r['node']==name and r['status']=='RUNNING' for r in active.values())
    cpu=int(n['CPUTot'])-int(n['CPUAlloc'])+own*64
    mem=int(n['RealMemory'])-int(n['AllocMem'])+own*192*1024
    count=max(0,min(cpu//64,mem//(192*1024)))
    slots[group]+=count
    node_slots.append(dict(node=name, group=group, slots_including_own_running=count, other_cpu_allocation=int(n['CPUAlloc'])-own*64))
slots['DCC']=6

destination={
    '4u3s':'fennario','5dc4':'fennario','5it3':'fennario',
    '5d68':'compsci','2rl0':'compsci','1a0r':'compsci','3cal':'compsci',
    '4pxf':'compsci','5a6y':'compsci','4wyq':'compsci','1gwc':'compsci',
    '2xgy':'DCC','3ma2':'DCC','4wyu':'DCC','2rfd':'DCC','3bua':'DCC','4wwi':'DCC'}
rows=[]
for r in new_designs:
    rows.append(dict(design_id=r['design_id'], system=r['system'], kind='new_tier', destination=destination[r['system']],
                     requested_cpus=64, limit_hours=336, state='PROPOSED_NOT_SUBMITTED'))
for r in inventory['remaining']:
    rows.append(dict(design_id=r['design'], system=r['system'], kind='historical_14day_noncompletion', destination=destination[r['system']],
                     requested_cpus=64, limit_hours=336, state='PROPOSED_NOT_SUBMITTED'))
assert len(rows)==len({r['design_id'] for r in rows})==24
assert not {r['design_id'] for r in rows} & (existing_local|existing_dcc)
assert all(r['destination']==active[d]['group'] for r in rows for d in active if d.split('_flex_')[0]==r['system'])
# A balanced 5/11/8 plan requires transferring the sole lower 4wwi tier too.
# This is only a proposal. Do not hold, cancel, or resubmit that pending case.
carry='4wwi_flex_p1'
assert carry in existing_local and carry not in active
carry_job=tasks[carry]['job']
assert after[carry_job][0]=='PENDING' and after[carry_job][2]=='Resources'
counts=dict(collections.Counter(r['destination'] for r in rows))
assert counts=={'compsci':11,'fennario':5,'DCC':8},counts
budget_per_slot={g:counts[g]*336/slots[g] for g in counts}
with (here/'markstar24_proposed.tsv').open('w') as f:
    writer=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
strict_counts=dict(counts)
strict_counts['compsci']+=2
strict_counts['DCC']-=2
summary=dict(snapshot_epoch=time.time(), status='PLAN_ONLY', new_tiers=7, previous_noncompletion=17, total=24,
             slots=dict(slots), node_slots=node_slots, proposed_counts=counts,
             added_budget_hours_per_slot=budget_per_slot,
             interpretation='Budget arithmetic using 14-day caps, not runtime forecasts or calendar finish estimates.',
             strict_existing_assignment_counts=strict_counts,
             proposed_existing_work_transfer=dict(design=carry, job=carry_job, from_group='compsci', to_group='DCC',
                                                   current_status='PENDING_UNCHANGED', applied=False,
                                                   reason='Keep all 4wwi tiers on one DCC CPU model; lower tier has not started.'),
             constraints=['Do not move any started design.',
                          'Use the same actual CPU model for all newly launched tiers of a system.',
                          'DCC concurrency is six total across original, transferred, and new workloads.',
                          'DCC must provide a partition allowing 14-day tasks and enough slots of one CPU model per system.',
                          '12825151/12825152 are original 48-CPU filtered runs and do not replace the frontier timeout reruns.',
                          'Wait for PACK* full workload measurements before launching the seven new MARK* designs.',
                          'Keep existing short work progressing; initially use at most 2/4/3 slots for added long tasks on fennario/compsci/DCC.',
                          'Raise long-task concurrency only as original backlog clears; no independent six-worker DCC array.'],
             external_snapshot=dict(completed=24, running=6, accepted_new_transfers=9,
                                    source='User report 2026-10-06, no new completions', individual_cpu_models='unknown'),
             retired_duplicate_jobs=ids, running_jobs_preserved=sorted(running_before), audit_root=str(out))
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
(here/'markstar_plan.json').write_text(json.dumps(summary, indent=2)+'\n')
registry=json.loads((here/'registry.json').read_text())
gpu_job=subprocess.check_output(['scontrol','show','job','12826506','-o'],text=True)
(out/'packstar_job.txt').write_text(gpu_job)
registry.update(packstar_job='12826506', packstar_status='SUBMITTED', packstar_limit_hours=6,
                packstar_grouped=True, packstar_cpus=128, packstar_gpus=4,
                packstar_output='/usr/xtmp/lz280/frontier_extension7_20261006/run_12826506',
                markstar_plan=str(here/'markstar_plan.json'))
(here/'registry.json').write_text(json.dumps(registry,indent=2)+'\n')
print('RETIRED_LOCAL_DUPLICATES', ids, flush=True)
print('SLOTS', dict(slots), flush=True)
print('PROPOSED_COUNTS',counts,'STRICT_COUNTS',strict_counts,flush=True)
print('PROPOSED_EXTRA_TRANSFER',carry,carry_job,'NOT_APPLIED',flush=True)
for g in ['fennario','compsci','DCC']:
    print(g, [r['design_id'] for r in rows if r['destination']==g],flush=True)
print('OUTPUT',out,flush=True)
