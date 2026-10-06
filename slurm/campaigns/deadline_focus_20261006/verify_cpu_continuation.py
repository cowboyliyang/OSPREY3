"""Verify live CPU dependency chains, frozen ungated runner, and DCC counts."""
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here=Path(__file__).parent
def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))
edges={r['job']:r['parent'] for r in read(here/'local_lane_dependencies.tsv')}
assert len(edges)==18 and len(set(edges.values()))==18, 'Each predecessor has exactly one successor'
roots={'12812439_'+str(i):('fennario' if i<4 else 'compsci') for i in [0,1,2,3,4,5,7,8,9,11,12,13]}
counts=collections.Counter()
for job in edges:
    cursor=job;seen=set()
    while cursor in edges:
        assert cursor not in seen,'Dependency cycle'
        seen.add(cursor);cursor=edges[cursor]
    assert cursor in roots,(job,cursor)
    counts[roots[cursor]]+=1
assert counts=={'compsci':13,'fennario':5}
raw=subprocess.check_output(['squeue','-r','-u','lz280','-h','-o','%i|%T|%E|%P|%C'],text=True)
queue={p[0]:p[1:] for line in raw.splitlines() if (p:=line.split('|'))}
for job,parent in edges.items():
    if job not in queue:continue
    state,dependency,partition,cpus=queue[job]
    assert cpus=='64' and '12826950' not in dependency,(job,queue[job])
    if state=='PENDING':assert dependency.startswith('afterany:'+parent), (job,dependency)
registry=json.loads((here/'registry.json').read_text())
frozen=Path(registry['markstar_prep'])/'markstar_runner.py'
content=frozen.read_text()
assert 'FOCUS_PACK_ROOT' not in content and 'screening.json' not in content
assert hashlib.sha256(frozen.read_bytes()).digest()==hashlib.sha256((here/'markstar_runner.py').read_bytes()).digest()
compile(content,str(frozen),'exec')
for path in here.glob('*.py'):compile(path.read_text(),str(path),'exec')
for path in here.glob('*.slurm'):subprocess.run(['bash','-n',str(path)],check=True)
remote=read(here/'dcc_actions.tsv')
assert all(r['gate']=='NONE' for r in remote)
remote_counts=dict(collections.Counter(r['action'] for r in remote))
assert remote_counts['DEFER_UNSTARTED']==57
assert remote_counts['KEEP_COVERAGE']+remote_counts['KEEP_PRIMARY']+remote_counts['KEEP_CONTINUATION_RESERVE']==10
priority=read(here/'dcc_launch_priority.tsv')
assert [int(r['priority']) for r in priority]==list(range(1,15))
future={r['design_id'] for r in remote if r['action'] in ('KEEP_COVERAGE','KEEP_PRIMARY','KEEP_CONTINUATION_RESERVE','PROPOSE_NEW_PRIMARY','PROPOSE_NEW_COVERAGE')}
assert len(priority)==len(future)==14 and {r['design_id'] for r in priority}==future
assert all(r['role'].startswith('LONG_') and r['gpu_dependency']=='NONE' for r in priority[:8])
assert all(r['role']=='COVERAGE_FALLBACK' for r in priority[8:])
retained={r['design_id'] for r in remote if r['action'] in ('KEEP_COVERAGE','KEEP_PRIMARY','KEEP_CONTINUATION_RESERVE')}
(here/'dcc_keep_existing_pending_cases.txt').write_text(''.join(r['design_id']+'\n' for r in priority if r['design_id'] in retained))
policy=json.loads((here/'policy.json').read_text())
policy['dcc_counts']=remote_counts
policy['dcc_launch_strategy']='LONG_FIRST_KEEP_SIX_CPU_JOBS_NO_GPU_WAIT'
policy['dcc_launch_priority']=[r['design_id'] for r in priority]
policy['dcc_long_candidates']=8
policy['dcc_short_or_unknown_coverage_fallbacks']=6
(here/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
subprocess.run(['git','diff','--check'],cwd=here.parents[2],check=True)
summary=dict(verification_job=os.environ['SLURM_JOB_ID'],local_future_designs=len(edges),by_cluster=dict(counts),
    lane_caps=dict(collections.Counter(roots.values())),all_active_future_jobs_request_64_cpus=True,
    frozen_runner_has_no_gpu_gate=True,dcc_gpu_gates_removed=True,dcc_existing_pending_keep=10,
    dcc_cancel_candidates=57,dcc_new=4,dcc_long_candidates_first=8,dcc_fallbacks_after_long_queue=6,
    current_running_roots=[j for j in roots if queue.get(j,[''])[0]=='RUNNING'])
(here/'cpu_continuation_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
