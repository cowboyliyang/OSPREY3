"""Record measured PACK* times; gate only the two newly proposed local targets."""
import csv
import json
import os
from pathlib import Path
import re

assert os.environ.get('SLURM_JOB_ID')
root=Path(os.environ['FOCUS_PACK_ROOT'])
records=list(csv.DictReader((root/'status.tsv').open(),delimiter='\t')) if (root/'status.tsv').exists() else []
result={}
for r in records:
    wall=Path(r['output'])/'wall.time' if r['output'] else None
    content=wall.read_text() if wall and wall.exists() else ''
    match=re.search(r'elapsed=(\d+(?:\.\d+)?)',content)
    seconds=float(match[1]) if match else None
    completed=r['status'] in ('COMPLETED','INCOMPLETE_ESTIMATES') and seconds is not None
    result[r['design']]=dict(status=r['status'],seconds=seconds,result_rows=r['result_rows'],
        estimated_rows=r['estimated_rows'],expected_rows=r['expected_rows'],
        eligible_for_primary_markstar=completed and seconds<=1800,output=r['output'])
summary=dict(pack_job=os.environ['FOCUS_PACK_JOB'],threshold_seconds=1800,
             note='Normal workload completion required. Returned binding estimates are recorded separately; only matched sequence results support numerical speedups.',designs=result)
(root/'screening.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
