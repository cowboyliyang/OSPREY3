"""Screen focused PACK* designs in one uninterrupted Slurm allocation."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

assert os.environ.get('SLURM_JOB_ID')
assert int(os.environ['SLURM_CPUS_PER_TASK']) == 128
prep = Path(os.environ['EXTENSION_PREP'])
assert (prep/'READY').is_file()
root = Path(os.environ['RUN_ROOT'])
package = prep/'package'
protocol = json.loads((package/'protocol.json').read_text())
started = time.time()
deadline = started + protocol['global_limit_seconds'] - protocol['shutdown_reserve_seconds']
with (root/'input_verification.log').open('w') as f:
    subprocess.run(['sha256sum', '-c', 'SHA256SUMS'], cwd=package, stdout=f, stderr=subprocess.STDOUT, check=True)
for name in ['portable_run.py', 'verification.json']:
    shutil.copy2(prep/name, root/name)
shutil.copy2(package/'designs.tsv', root/'designs.tsv')
(root/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
rows = list(csv.DictReader((root/'designs.tsv').open(), delimiter='\t'))
assert len(rows)==10 and len({r['design_id'] for r in rows})==10
selected = {int(i) for i in os.environ['FOCUS_PACK_INDICES'].split(',')}
assert selected and selected <= {int(r['task_id']) for r in rows}
rows = [r for r in rows if int(r['task_id']) in selected]
env = dict(os.environ, BUILD_ROOT=str(prep/'build'), INPUT_ROOT=str(package), RESULT_ROOT=str(root/'results'))
common = ['python3', str(root/'portable_run.py'), '--frontier-manifest', str(root/'designs.tsv'),
          '--arm', 'pair-only', '--rb', '1', '--seed', '42', '--gpus', '4', '--gpu-gib', '85']
records = []
fields = ['index','design','start','end','budget_seconds','runner_exit_code','status','result_rows','expected_rows','estimated_rows','output']
with (root/'status.tsv').open('w') as status_file:
    writer = csv.DictWriter(status_file, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    status_file.flush()
    for position, row in enumerate(rows):
        index, design = row['task_id'], row['design_id']
        case_start = time.time()
        reserve = (len(rows)-position-1)*protocol['future_case_reserve_seconds']
        budget = min(protocol['per_design_max_seconds'], int(deadline-case_start-reserve))
        info, output_path, code = {}, '', None
        status = 'NOT_STARTED_TIME_BUDGET'
        print('DESIGN_START', design, 'budget_seconds', budget, flush=True)
        if budget >= 60:
            with (root/('design_'+index+'.out')).open('w') as f:
                preflight = common+['--index', index, '--mode', 'preflight', '--heap-gib', '8', '--host-gib', '4']
                code = subprocess.run(['timeout', '--signal=TERM', '--kill-after=10s', str(min(300,budget))+'s', *preflight],
                                      env=env, stdout=f, stderr=subprocess.STDOUT).returncode
                status = 'PREFLIGHT_TIMEOUT' if code in (124,137) else 'PREFLIGHT_FAILED'
                if code==0:
                    remaining = max(1, int(min(case_start+budget, deadline)-time.time()))
                    command = common+['--index', index, '--mode', 'full', '--heap-gib', '850', '--host-gib', '800']
                    code = subprocess.run(['timeout', '--signal=TERM', '--kill-after=10s', str(remaining)+'s', *command],
                                          env=env, stdout=f, stderr=subprocess.STDOUT).returncode
                    manifests = list((root/'results').glob(design+'_full_*/manifest.json'))
                    if manifests:
                        assert len(manifests)==1
                        info = json.loads(manifests[0].read_text())
                        output_path = str(manifests[0].parent)
                    if code in (124,137):
                        status = 'TIMEOUT'
                    elif (info.get('exit_code')==0 and info.get('result_rows')==int(row['expected_sequences'])
                          and info.get('status') in ('COMPLETED','INCOMPLETE_ESTIMATES')
                          and (Path(output_path)/'wall.time').is_file()):
                        status = info['status']
                    else:
                        status = 'FAILED_OR_INCOMPLETE'
        record = dict(index=index, design=design, start=case_start, end=time.time(), budget_seconds=budget,
                      runner_exit_code=code, status=status, result_rows=info.get('result_rows',''),
                      expected_rows=row['expected_sequences'], estimated_rows=info.get('estimated_rows',''), output=output_path)
        writer.writerow(record)
        status_file.flush()
        records.append(record)
        (root/'summary.json').write_text(json.dumps(dict(job=os.environ['SLURM_JOB_ID'], prep=str(prep),
                                                       allocation_kept_between_designs=True, attempts=records), indent=2)+'\n')
        print('DESIGN_END', json.dumps(record), flush=True)
print('ALL_SELECTED_DESIGNS_ATTEMPTED', flush=True)
raise SystemExit(0 if all(r['status'] in ('COMPLETED','INCOMPLETE_ESTIMATES') for r in records) else 5)
