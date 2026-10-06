"""Append validated workloads without changing any existing Slurm job."""
import collections
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID')
here = Path(__file__).parent
previous = here.parent/'deadline_focus_20261006'
registry = json.loads((here/'registry.json').read_text())
assert not registry.get('submissions'), 'Inspect existing submissions; do not replay'
prep = Path(registry['prep_root'])
assert (prep/'READY').is_file()
out = prep.parent/('launch_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)

def read(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))

def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, restval='NA', delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def save_registry():
    (here/'registry.json').write_text(json.dumps(registry, indent=2)+'\n')

def submit(args, design):
    job = subprocess.check_output(['sbatch','--parsable','--account=grisman',*args], text=True).strip().split(';')[0]
    assert job.isdigit(), job
    record = dict(job=job, design=design, arguments=args)
    registry.setdefault('submissions', []).append(record)
    save_registry()
    with (out/'submitted_jobs.jsonl').open('a') as f:
        f.write(json.dumps(record)+'\n')
    return job

rows = read(prep/'package/designs.tsv')
assert len(rows) == 9 and all(r['pack_platform']=='gpu' for r in rows)
assert '2p4a_flex_p12' in {r['design_id'] for r in rows}
assert '2p4a_flex_p18' not in {r['design_id'] for r in rows}
assert read(here/'requested_designs.tsv') == read(prep/'package/provenance/requested_designs.tsv')
subprocess.run(['sha256sum','-c','SHA256SUMS'], cwd=prep/'package', check=True, stdout=subprocess.DEVNULL)
for path in here.glob('*.py'):
    compile(path.read_text(), str(path), 'exec')
for path in here.glob('*.slurm'):
    subprocess.run(['bash','-n',str(path)], check=True)

# Existing serial lanes stay intact, including their smaller coverage cases.
# Add at lane tails, not alongside existing successors, to preserve the cap of 8.
parents = {'1b6c_flex_p4':'12827032_2', '2hnu_flex_p11':'12827031_1',
           '3k3q_flex_p3':'12826674_9', '4znc_flex_p14':'12827030_0'}
existing_edges = read(previous/'local_lane_dependencies.tsv')
assert not set(parents.values()) & {r['parent'] for r in existing_edges}
assert len(set(parents.values())) == 4
queue = subprocess.check_output(['squeue','-u','lz280','-r','-h','-o','%i|%T|%R'], text=True)
(out/'queue_before.txt').write_text(queue)
for parent in parents.values():
    detail = subprocess.check_output(['scontrol','show','job',parent,'-o'], text=True)
    assert 'JobState=PENDING' in detail, (parent, detail)
    (out/(parent+'.before.txt')).write_text(detail)

paper = json.loads(Path('/home/users/lz280/BranchMARK*_paper/recomb_draft/runtime_records/data.json').read_text())['frontier']
historical = {r['design_id']:r for r in paper['runs'] if r['method']=='MARK*'}
markprep = out/'markstar'
package = markprep/'package'
package.mkdir(parents=True)
configs, markrows = {}, []
indices = collections.Counter()
for r in rows:
    source = Path(historical[r['anchor_design']]['path'])
    command = shlex.split((source/'command.sh').read_text())
    changes = {'-XX:ActiveProcessorCount=':'64', '-Dosprey.bench.numCPUs=':'64',
               '-Dosprey.bench.designId=':r['design_id'], '-Dosprey.bench.outputDir=':'{RUN}',
               '-Djava.io.tmpdir=':'{RUN}/tmp', '-Dosprey.bench.pdbPath=':str(prep/'package'/r['pdb_relative']),
               '-Dosprey.bench.mutable=':r['mutable'], '-Dosprey.bench.flexible=':r['flexible']}
    rewritten, changed = [], set()
    for token in command:
        key = next((k for k in changes if token.startswith(k)), None)
        if key:
            changed.add(key)
            token = key+changes[key]
        rewritten.append(token.replace(str(source),'{RUN}'))
    assert changed == set(changes)
    assert '-Dosprey.bench.method=markstar' in rewritten
    assert all(Path(p).exists() for p in rewritten[rewritten.index('-cp')+1].split(':'))
    heap = int(re.fullmatch(r'-Xmx(\d+)g', next(t for t in rewritten if t.startswith('-Xmx')))[1])
    assert heap < 192
    group = r['destination']
    row = dict(r, group_index=indices[group], requested_cpus=64, requested_memory_gib=192,
               java_heap_gib=heap, limit_hours=336, reference_design=r['anchor_design'],
               reference_mark_path=str(source), source_package=str(prep/'package'))
    indices[group] += 1
    markrows.append(row)
    configs[r['design_id']] = dict(row=row, command_template=rewritten,
        source_command_sha256=hashlib.sha256((source/'command.sh').read_bytes()).hexdigest(),
        model_contains='AMD EPYC 9554' if group=='compsci' else None,
        input_pdb=str(prep/'package'/r['pdb_relative']))
write(package/'designs.tsv', markrows)
(package/'run_configs.json').write_text(json.dumps(configs, indent=2)+'\n')
for name in ['markstar_runner.py','markstar_array.slurm']:
    shutil.copy2(here/name, markprep/name)
assert 'FOCUS_PACK_ROOT' not in (markprep/'markstar_runner.py').read_text()
(markprep/'READY').write_text('Nine MARK* and GPU PACK* input preflights passed\n')
registry['markstar_prep'] = str(markprep)
registry['launch_audit'] = str(out)
save_registry()

# The old twelve-design GPU allocation runs first. No CPU PACK* job is submitted.
pack = submit(['--dependency=afterany:12827169', '--export=ALL,EXTENSION_PREP='+str(prep),
               str(prep/'run_group.slurm')], 'ALL_NINE_GPU_PACKSTAR')
packroot = prep.parent/('pack_'+pack)
registry.update(packstar_job=pack, packstar_output=str(packroot), packstar_status='SUBMITTED',
                packstar_after_job='12827169', packstar_cpus=128, packstar_gpus=4, packstar_limit_hours=6)
save_registry()
screen = submit(['--dependency=afterany:'+pack, '--export=ALL,PACK_RUN_ROOT='+str(packroot),
                 str(here/'collect_pack.slurm')], 'PACKSTAR_MEASUREMENT_AUDIT')
registry['packstar_screening_job'] = screen
save_registry()

local_jobs = []
for row in markrows:
    if row['destination'] != 'compsci':
        continue
    index = str(row['group_index'])
    parent = parents[row['design_id']]
    job = submit(['--job-name=markstar-near14-compsci', '--partition=compsci', '--exclude=linux[31-40]',
                  '--array='+index+'-'+index+'%1', '--dependency=afterany:'+parent,
                  '--export=ALL,MARKSTAR_EXTENSION_PREP='+str(markprep)+',MARKSTAR_EXTENSION_GROUP=compsci',
                  str(markprep/'markstar_array.slurm')], row['design_id'])
    local_jobs.append(dict(design=row['design_id'], job=job+'_'+index, group='compsci',
                           after_job=parent, gpu_dependency='NONE', prep=str(markprep)))
    registry['local_markstar_jobs'] = local_jobs
    save_registry()
write(here/'local_submissions.tsv', local_jobs)
write(here/'dcc_additions.tsv', [r for r in markrows if r['destination']=='DCC'])
(here/'dcc_additions_cases.txt').write_text(''.join(r['design_id']+'\n' for r in markrows if r['destination']=='DCC'))

# Extend the active DCC handoff without changing its first fourteen priorities.
dcc_rows = read(previous/'dcc_new_designs.tsv')
assert len(dcc_rows)==4 and not {r['design_id'] for r in dcc_rows} & {r['design_id'] for r in markrows}
dcc_rows += [r for r in markrows if r['destination']=='DCC']
for index, row in enumerate(dcc_rows):
    row['handoff_index'] = index
write(previous/'dcc_new_designs.tsv', dcc_rows)
(previous/'dcc_new_cases.txt').write_text(''.join(r['design_id']+'\n' for r in dcc_rows))
priorities = read(previous/'dcc_launch_priority.tsv')
assert len(priorities)==14
for row in markrows:
    if row['destination']=='DCC':
        priorities.append(dict(priority=len(priorities)+1, design_id=row['design_id'],
            role='CONSERVATIVE_ADDITION' if row['system']=='2p4a' else 'NEAR14_EXPLORATORY_ADDITION',
            prior_submission='NEW_DEFINITION_NOT_SUBMITTED_REMOTELY', gpu_dependency='NONE',
            runtime_evidence='Unmeasured; preserve all prior smaller cases; no guaranteed near-14-day runtime'))
write(previous/'dcc_launch_priority.tsv', priorities)
dcc_actions = read(previous/'dcc_actions.tsv')
for r in dcc_actions:
    if r['design_id']=='2xxm_flex_p1':
        r.update(action='KEEP_COMPLETED', reported_status='COMPLETED_USER_REPORT')
for row in markrows:
    if row['destination']=='DCC':
        dcc_actions.append(dict(design_id=row['design_id'], action='APPEND_EXPLORATORY_ADDITION',
            reported_status='HANDOFF_NOT_SUBMITTED_REMOTELY', requested_cpus=64, limit_hours=336,
            combined_concurrency=6, gate='NONE',
            cpu_rule='Verify and use the actual CPU model of the current-round same-system anchor; do not infer from partition'))
write(previous/'dcc_actions.tsv', dcc_actions)

external = json.loads((previous/'external_status_snapshot.json').read_text())
external['completed_seconds']['2xxm_flex_p1'] = 6338
external['reported_completed_count'] = len(external['completed_seconds'])
external['completion_correction'] = json.loads((here/'selection.json').read_text())['external_correction']
(previous/'external_status_snapshot.json').write_text(json.dumps(external, indent=2)+'\n')
actions = read(previous/'design_actions.tsv')
for r in actions:
    if r['design']=='2xxm_flex_p1':
        r.update(action='KEEP_COMPLETED', current_status='COMPLETED_DCC_USER_REPORT')
local_by_name = {r['design']:r['job'] for r in local_jobs}
for row in markrows:
    r = {key:'NA' for key in actions[0]}
    r.update(design=row['design_id'], system=row['system'], action='APPEND_EXPLORATORY_ADDITION',
             current_local_job=local_by_name.get(row['design_id'], ''),
             current_status='SUBMITTED_LOCAL' if row['destination']=='compsci' else 'HANDOFF_NOT_SUBMITTED_REMOTELY')
    # Preserve the active inventory's column names and original rows.
    for key in ('current_assignment','proposed_destination'):
        if key in r:
            r[key] = 'local' if row['destination']=='compsci' else 'dcc'
    actions.append(r)
write(previous/'design_actions.tsv', actions)
coverage = read(previous/'system_coverage.tsv')
for r in coverage:
    if r['system']=='2xxm':
        r.update(coverage_state='COMPLETED', completed='2xxm_flex_p1', pending_or_new='')
    new = next((row for row in rows if row['system']==r['system']), None)
    if new:
        r['pending_or_new'] = ';'.join(filter(None, [r['pending_or_new'],new['design_id']]))
        r['primary'] = new['design_id']
write(previous/'system_coverage.tsv', coverage)
edges = existing_edges + [dict(job=r['job'], parent=r['after_job'], action='APPEND_AFTER_EXISTING_LANE') for r in local_jobs]
write(previous/'local_lane_dependencies.tsv', edges)
policy = json.loads((previous/'policy.json').read_text())
policy['coverage_unconfirmed'].pop('2xxm_flex_p1', None)
policy['dcc_counts'] = dict(collections.Counter(r['action'] for r in dcc_actions))
policy['dcc_launch_priority'] = [r['design_id'] for r in priorities]
policy['near14_additions'] = dict(registry=str(here/'registry.json'), mode='APPEND_ONLY_SMALLER_FIRST',
    designs=[r['design_id'] for r in rows], all_packstar_on_gpu=True, packstar_job=pack,
    local_markstar_jobs=local_jobs, dcc_status='HANDOFF_ONLY_NOT_APPLIED_REMOTELY',
    two_xxm_p1_completed=True, conservative_target='2p4a_flex_p12', withdrawn_target='2p4a_flex_p18')
policy['rules'].append('October 6 final revision: preserve the selected smaller queue; append nine designs, use GPU PACK* for all nine, and use conservative 2p4a p12 rather than p18.')
(previous/'policy.json').write_text(json.dumps(policy, indent=2)+'\n')
old_registry = json.loads((previous/'registry.json').read_text())
old_registry['near14_additions_registry'] = str(here/'registry.json')
(previous/'registry.json').write_text(json.dumps(old_registry, indent=2)+'\n')
registry.update(markstar_status='FOUR_LOCAL_SUBMITTED_FIVE_DCC_HANDOFF_ONLY', dcc_status='NOT_APPLIED_REMOTELY',
    dcc_combined_concurrency=6, local_lane_limit=8, existing_jobs_modified=[],
    applied_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save_registry()
(out/'summary.json').write_text(json.dumps(registry, indent=2)+'\n')
print('SUBMITTED', json.dumps(dict(packstar_job=pack, local_markstar_jobs=local_jobs,
                                 dcc_additions=5, dcc_status=registry['dcc_status'])), flush=True)
