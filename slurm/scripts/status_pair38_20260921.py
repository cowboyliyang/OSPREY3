"""Read-only run census; write a separate interim report through Slurm."""
import os,json,csv,collections,math
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
root=Path('/usr/xtmp/lz280/packstar_pair38_20260921/launch_12677939')
source=(root/'source/pair38_20260921.py').read_text()
old="report = root / 'comparison'"
assert source.count(old)==1
source=source.replace(old,"report = root / ('interim_' + os.environ['SLURM_JOB_ID'])")
task_read="task = json.loads(task_file.read_text()) if task_file.exists() else {}"
assert source.count(task_read)==1
source=source.replace(task_read,task_read+"\n        task = TASK_OVERRIDES.get(case['index'], task)")
overrides={};audits=[]
namespace={'__name__':'pair38_interim','__file__':str(root/'source/pair38_20260921.py'),'TASK_OVERRIDES':overrides}
exec(compile(source,'pair38_interim_snapshot.py','exec'),namespace)
# The original audit incorrectly allowed only Estimated/Aborted. Unstable is
# a legitimate non-success terminal status. Reaudit independently; frozen task
# records and calculation results are preserved verbatim.
plan=json.loads((root/'plan.json').read_text())
for case in plan['cases']:
    path=root/'tasks'/(str(case['index'])+'.json')
    if not path.exists():continue
    task=json.loads(path.read_text())
    if task.get('status')!='AUDIT_FAILED':continue
    if task['design']!='3bu8' or "'prot_status': 'Unstable'" not in task.get('error',''):continue
    dest=Path(task['run']);manifest=json.loads((dest/'manifest.json').read_text())
    prior=json.loads((root/'baseline'/case['design']/'manifest.json').read_text())
    assert task['java_exit']==manifest['exit_code']==0 and task['runner_exit'] in (0,5)
    assert manifest['pdb_sha256']==prior['pdb_sha256']==case['pdb_sha256']
    props=manifest['properties'];old_props=prior['properties']
    different={k for k in set(props)|set(old_props) if props.get(k)!=old_props.get(k)}
    assert not different-namespace['ALLOWED_PROPERTY_CHANGES']
    assert props['branchdp.cutoff.residualBudget']=='1' and props['packstar.pac.frequencySeverity.tripleEta']=='false'
    assert props['packstar.pac.frequencySeverity.jointMomentLearning']=='true'
    rr=namespace['rows'](dest/(case['design']+'_packstar.csv'))
    assert len(rr)==case['expected_rows'] and {r['sequence'].strip() for r in rr}==set(case['expected_sequences'])
    protocols=0
    for state in (dest/'adaptive_frequency_severity').glob('state-*'):
        protocol=namespace['keys'](state/'frequency_severity_protocol.tsv')
        if protocol:
            assert protocol['tripleEtaEnabled']=='false' and protocol['proposalLearningEnabled']=='true'
            protocols+=1
        for p in list(state.glob('*selection.tsv'))+[state/'frequency_severity_final_interval.tsv']:
            data=namespace['keys'](p)
            assert int(data.get('selectedTripleEtaPositionTriples','0'))==0
    assert protocols>0
    states=collections.Counter()
    for r in rr:
        for state in ('prot','lig','comp'):
            status=r[state+'_status'];states[status]+=1
            assert status in ('Estimated','Aborted','Unstable')
            if status=='Estimated':assert math.isfinite(float(r[state+'_eps'])) and float(r[state+'_eps'])<=.683001
    overrides[case['index']]=dict(task,status='COMPLETE_REAUDITED',technical_ok=True)
    audits.append(dict(design=case['design'],original_error=task['error'],state_statuses=dict(states),
        decision='Include completed calculation; Unstable remains non-success, preserve raw state statuses',protocols=protocols))
namespace['summarize'](root)
report=root/('interim_'+os.environ['SLURM_JOB_ID'])
(report/'analysis_source.py').write_text(Path(__file__).read_text())
(report/'summarizer_snapshot.py').write_text(source)
(report/'independent_reaudit.json').write_text(json.dumps(audits,indent=2)+'\n')
def rows(name):
    with (report/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
systems=rows('systems.tsv');seqs=rows('sequences.tsv')
complete={r['design'] for r in systems if r['technical_ok']=='True'}
observed=[r for r in seqs if r['design'] in complete]
matched=[r for r in observed if r['triple_baseline'] in ['Estimated','Aborted'] and r['pair_only'] in ['Estimated','Aborted']]
result=dict(completed_systems=len(complete),remaining_systems=[r['design'] for r in systems if r['design'] not in complete],
    observed_sequences=len(observed),pair_success_observed=sum(r['pair_only']=='Estimated' for r in observed),
    matched_sequences=len(matched),matched_pair_success=sum(r['pair_only']=='Estimated' for r in matched),
    matched_triple_success=sum(r['triple_baseline']=='Estimated' for r in matched),
    matched_transitions=dict(collections.Counter(r['triple_baseline']+' -> '+r['pair_only'] for r in matched)),
    pair_all_sequences_success_systems=[r['design'] for r in systems if r['technical_ok']=='True' and int(r['pair_estimated'])==int(r['expected'])])
for k,denom in [('pair_success_observed','observed_sequences'),('matched_pair_success','matched_sequences'),('matched_triple_success','matched_sequences')]:
    result[k+'_percent']=100*result[k]/result[denom] if result[denom] else None
(report/'status.json').write_text(json.dumps(result,indent=2)+'\n')
print('INTERIM_STATUS',json.dumps(result,indent=2));print('REPORT',report)
