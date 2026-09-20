"""Write a small versioned design-to-job registry after authorized submission."""
import csv, hashlib, json, os, sys
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
prep=Path(sys.argv[1]); jobs=dict(zip(['compsci','grisman','fennario'],sys.argv[2:5]))
assert len(jobs)==3 and all(x.isdigit() for x in jobs.values())
assert (prep/'READY').read_text().strip()=='52 designs validated'
repo=Path('/home/users/lz280/OSPREY3-fresh-packstar')
for line in (prep/'package/SHA256SUMS').read_text().splitlines():
    checksum,relative=line.split('  ',1)
    assert hashlib.sha256((prep/'package'/relative).read_bytes()).hexdigest()==checksum
designs=json.loads((prep/'package/designs.json').read_text())
fields=next(csv.reader((repo/'slurm/h200/frontier.tsv').open(),delimiter='\t'))
rows=[]
for r in designs:
    row={k:r.get(k,'') for k in fields}
    row.update(markstar_job=f'{jobs[r["scheduling_group"]]}_{r["task_id"]}',markstar_cpus=16,source_package=str(prep/'package'),source_task_id=r['task_id'])
    rows.append(row)
target=repo/'slurm/h200/frontier_remaining18_plus4.tsv'
assert not target.exists()
with target.open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
audit=dict(prep_root=str(prep),jobs=jobs,designs=52,systems=22,cpus=16,memory_gib=96,heap_gib=64,time_limit_days=14,account='grisman',normalizations=json.loads((prep/'package/config/input_normalizations.json').read_text()))
(repo/'slurm/h200/frontier_remaining18_plus4.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit),flush=True)
