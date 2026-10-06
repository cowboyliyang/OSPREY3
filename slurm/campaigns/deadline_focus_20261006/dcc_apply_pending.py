"""Apply the exact DCC defer list to a reviewed design_id/job_id TSV.

Run inside a Slurm management job with --account=grisman. This never submits
benchmarks or edits an existing array's manifest. Default is a dry run.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Submit this management script through Slurm'
p=argparse.ArgumentParser()
p.add_argument('--job-map',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--apply',action='store_true')
args=p.parse_args()
here=Path(__file__).parent
cancel=set((here/'dcc_cancel_if_pending_cases.txt').read_text().splitlines())
protected=set((here/'dcc_protected_cases.txt').read_text().splitlines())
assert not cancel & protected
mapping=list(csv.DictReader(args.job_map.open(),delimiter='\t'))
assert len({r['job_id'] for r in mapping})==len(mapping), 'Job IDs must identify individual tasks'
args.output.mkdir(parents=True,exist_ok=False)
records=[]
for r in mapping:
    design,job=r['design_id'],r['job_id']
    assert re.fullmatch(r'\d+(?:_\d+)?',job),job
    if design not in cancel:continue
    query=subprocess.run(['scontrol','show','job',job,'-o'],text=True,capture_output=True)
    (args.output/(job+'.txt')).write_text(query.stdout+query.stderr)
    info=dict(re.findall(r'(\w+)=([^\s]+)',query.stdout))
    state=info.get('JobState','NOT_IN_LIVE_QUEUE')
    action='SKIP_STARTED_OR_FINISHED'
    # Array parents/ranges are prohibited, even when only one row maps to them.
    individual=('ArrayJobId' not in info or info.get('ArrayTaskId','').isdigit() and '_' in job)
    never_started=info.get('Restarts')=='0' and info.get('RunTime')=='00:00:00'
    if state=='PENDING' and individual and never_started:
        action='CANCEL_PENDING' if args.apply else 'WOULD_CANCEL_PENDING'
        if args.apply:subprocess.run(['scancel','--state=PENDING',job],check=True)
    records.append(dict(design_id=design,job_id=job,state_before=state,action=action))
    (args.output/'actions.json').write_text(json.dumps(records,indent=2)+'\n')
missing=sorted(cancel-{r['design_id'] for r in mapping})
(args.output/'unmapped_designs.json').write_text(json.dumps(missing,indent=2)+'\n')
print(json.dumps(dict(apply=args.apply,actions=records,unmapped_designs=missing),indent=2))
