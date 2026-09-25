"""Frozen pair-only RB=2/4 launch, using the RB=1 protocol and Slurm audits."""
import argparse
import ast
from collections import Counter
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pair38_rb24_kernel_20260922 as kernel

REPO = kernel.REPO
PARENT = Path('/usr/xtmp/lz280/packstar_pair38_rb24_20260922')
RB1 = Path('/usr/xtmp/lz280/packstar_pair38_20260921/launch_12677939')
HISTORY = Path('/usr/xtmp/lz280/packstar_pair38_rb24_20260922/audit_12681033/history.json')
BASE = Path('/usr/xtmp/lz280/packstar_rb38_20260918/launch_12637884')
NAMES = ('pair38_rb24_20260922.py', 'pair38_rb24_kernel_20260922.py',
         'launch_pair38_rb24_20260922.slurm', 'run_pair38_rb24_20260922.slurm',
         'summarize_pair38_rb24_20260922.slurm')


def scheduler(root):
    text = subprocess.check_output(['scontrol', 'show', 'nodes', '--oneliner'], text=True)
    (root / 'scheduler_before_submission.txt').write_text(text)
    nodes = []
    for line in text.splitlines():
        d = dict(re.findall(r'(\w+)=([^ ]*)', line))
        if not set(d.get('Partitions', '').split(',')).intersection({'grisman','compsci-gpu'}):
            continue
        if any(s in d.get('State','') for s in ('DOWN','DRAIN','NOT_RESPONDING')):
            continue
        def tres(value):
            return dict(p.split('=',1) for p in value.split(',') if '=' in p)
        cfg, alloc = tres(d.get('CfgTRES','')), tres(d.get('AllocTRES',''))
        free = {k[9:]:int(v)-int(alloc.get(k,0)) for k,v in cfg.items() if k.startswith('gres/gpu:')}
        nodes.append(dict(node=d['NodeName'], partition=d['Partitions'],
            cpus=int(d['CPUTot'])-int(d['CPUAlloc']),
            mem_gib=(int(d['RealMemory'])-int(d['AllocMem']))/1024,
            gpus=free))
    kernel.dump(root / 'available_resources.json', nodes)
    return nodes


def allocate(design, history, old, nodes):
    # Same requested CPU/GPU counts and host budgets for both new RB values.
    seconds = max(h['elapsed_s'] for h in history if h['rb'] in (2,4))
    if design in ('3bua','4z80'):
        tier, cpus, gpus = 'very_hard', 64, 8
    elif seconds >= 15000:
        tier, cpus, gpus = 'hard', 48, 4
    elif seconds >= 7000:
        tier, cpus, gpus = 'medium', 32, 2
    elif seconds >= 2000:
        tier, cpus, gpus = 'light', 16, 1
    else:
        tier, cpus, gpus = 'easy', 8, 1
    resources = {k:old[k] for k in ('mem_gib','heap_gib','host_gib')}
    # Sparse 2xxm completed in under 1.5h with <7GiB RSS. Keep 144GiB
    # host headroom, but avoid inheriting RB=1's 512GiB reservation.
    if design == '2xxm':
        resources = dict(mem_gib=192, heap_gib=160, host_gib=144)
    # Use GPU families already exercised by the frozen RB=1 build. Leave
    # GPU placement to Slurm so an intervening allocation cannot strand jobs.
    families = ['a5000'] if gpus >= 4 else ['titan_v','rtx_5000','p100','v100','a5000']
    candidates = []
    for node in nodes:
        if node['cpus'] < cpus or node['mem_gib'] < resources['mem_gib']:
            continue
        for family in families:
            if node['gpus'].get(family,0) >= gpus:
                candidates.append((families.index(family), node['cpus']-cpus, node['node'], family, node))
    if candidates:
        _, _, candidate, family, node = min(candidates, key=lambda item:item[:4])
        node['cpus'] -= cpus
        node['mem_gib'] -= resources['mem_gib']
        node['gpus'][family] -= gpus
        partition = 'grisman' if family == 'titan_v' or gpus == 8 else (
            'grisman,compsci-gpu' if family == 'a5000' else 'compsci-gpu')
        availability = 'fits_snapshot'
    else:
        candidate, family = 'scheduler', 'a5000'
        partition = 'grisman' if gpus == 8 else 'grisman,compsci-gpu'
        availability = 'queued_capacity'
    return dict(resources, tier=tier, cpus=cpus, gpus=gpus, gpu_type=family,
        partition=partition, node='scheduler', snapshot_candidate=candidate,
        snapshot_availability=availability, difficulty_seconds=seconds,
        time_limit='7-00:00:00' if design in ('3bua','4z80') else '3-00:00:00')


def launch():
    root = PARENT / ('launch_' + os.environ['SLURM_JOB_ID'])
    source = root / 'source'
    source.mkdir(parents=True, exist_ok=False)
    for name in NAMES:
        shutil.copy2(REPO / 'slurm/scripts' / name, source / name)
    shutil.copy2(RB1 / 'source/portable_run.py', source / 'portable_run.py')
    shutil.copy2(HISTORY, root / 'historical_workloads.json')
    for path in source.glob('*.py'):
        ast.parse(path.read_text())
    for path in source.glob('*.slurm'):
        subprocess.run(['bash','-n',str(path)], check=True)
    assert (kernel.BUILD / 'READY').is_file()
    for i, relative in enumerate(('src/main', 'src/test/java/edu/duke/cs/osprey/markstar/bench/GenericPDBBench.java',
                                  'slurm/h200/production.properties', 'slurm/h200/baseline38.csv')):
        kernel.checked(['diff','-qr',str(kernel.BASE_BUILD / 'source' / relative),str(kernel.BUILD / 'source' / relative)],
            root / f'baseline_source_identity_{i}.txt')
        # The working tree changed after RB1. Record it, but deliberately
        # execute the exact frozen RB1 build, never the unbuilt working tree.
        with (root / f'current_source_difference_{i}.txt').open('w') as stream:
            status = subprocess.run(['diff','-qr',str(REPO/relative),str(kernel.BUILD/'source'/relative)],
                stdout=stream,stderr=subprocess.STDOUT)
            assert status.returncode in (0,1)
    reference = json.loads((RB1 / 'plan.json').read_text())
    history = json.loads(HISTORY.read_text())
    assert len(reference['cases']) == 38 and reference['total_sequences'] == 1103
    assert reference['build'] == str(kernel.BUILD)
    assert kernel.sha(source/'portable_run.py') == reference['source_sha256']['portable_run.py']
    free = shutil.disk_usage(PARENT).free
    assert free > 40 * 2**40, f'Insufficient xtmp capacity for conservative 38TiB temporary estimate: {free}'
    kernel.dump(root/'storage_plan.json', dict(free_bytes=free, download_bytes=0,
        input_source_allowance_bytes=50*2**20, output_allowance_bytes=200*2**30,
        temporary_upper_allowance_bytes=38*2**40, approximate_file_allowance=200000))
    nodes = scheduler(root)
    groups = {}
    all_cases = []
    for rb in (2,4):
        group = root / f'rb{rb}'
        for name in ('source','inputs','baseline','tasks','runs'):
            (group / name).mkdir(parents=True)
        for path in source.iterdir():
            shutil.copy2(path, group/'source'/path.name)
        cases = []
        for old in reference['cases']:
            design = old['design']
            h = next(h for h in history if h['rb']==rb and h['design']==design)
            prior = Path(h['path'])
            manifest = json.loads((prior/'manifest.json').read_text())
            ref_manifest = json.loads((RB1/'baseline'/design/'manifest.json').read_text())
            props, ref_props = manifest['properties'], ref_manifest['properties']
            different = {key for key in set(props)|set(ref_props) if props.get(key)!=ref_props.get(key)}
            assert not different - kernel.ALLOWED_PROPERTY_CHANGES - {'branchdp.cutoff.residualBudget'}, (design, different)
            assert props['branchdp.cutoff.residualBudget'] == str(rb)
            assert props['packstar.pac.frequencySeverity.jointMomentLearning'] == 'true'
            assert props['packstar.pac.frequencySeverity.tripleEta'] == 'true'
            assert manifest['options']['index'] == old['index'] and manifest['options']['seed'] == 42
            target = group/'inputs'/design/Path(old['pdb']).name
            target.parent.mkdir()
            shutil.copy2(old['pdb'], target)
            assert kernel.sha(target) == old['pdb_sha256'] == manifest['pdb_sha256']
            evidence = group/'baseline'/design
            evidence.mkdir()
            for name in ('manifest.json','wall.time','gpus.csv',design+'_packstar.csv'):
                if (prior/name).exists():
                    shutil.copy2(prior/name, evidence/name)
            prior_rows = kernel.rows(evidence/(design+'_packstar.csv'))
            if prior_rows:
                assert len(prior_rows) == old['expected_rows']
                assert {r['sequence'].strip() for r in prior_rows} == set(old['expected_sequences'])
            cases.append(dict(index=old['index'], design=design, rb=rb,
                pdb=str(target), pdb_sha256=old['pdb_sha256'], mutable=old['mutable'], flexible=old['flexible'],
                expected_rows=old['expected_rows'], expected_sequences=old['expected_sequences'],
                baseline=str(prior), historical_elapsed_s=h['elapsed_s'], historical_status=h['status'],
                historical_peak_rss_gib=h['peak_rss_gib'],
                baseline_technical_note='TIMEOUT at 2 days (sacct 12638022); no final CSV' if rb==4 and design=='4z80' else '',
                group_root=str(group)))
        groups[rb] = cases
        all_cases += cases
    old_by_design = {c['design']:c for c in reference['cases']}
    def difficulty(c):
        if c['design']=='4z80':
            return 172804
        return max(h['elapsed_s'] for h in history if h['design']==c['design'] and h['rb'] in (2,4))
    # Interleave RB values by system difficulty; neither arm waits for the
    # other arm's entire cohort. Slurm places both against shared capacity.
    all_cases.sort(key=lambda c:(-difficulty(c),c['design'],c['rb']))
    for case in all_cases:
        hh = [h for h in history if h['design']==case['design']]
        case.update(allocate(case['design'], hh, old_by_design[case['design']], nodes))
    assert len(all_cases) == len({(c['rb'],c['index']) for c in all_cases}) == 76
    addon_root = root/'rb4_triple_4z80'
    for name in ('source','inputs','baseline','tasks','runs'):
        (addon_root/name).mkdir(parents=True)
    for path in source.iterdir():
        shutil.copy2(path,addon_root/'source'/path.name)
    triple_case = dict(next(c for c in groups[4] if c['design']=='4z80'))
    triple_case['group_root'] = str(addon_root)
    pdb = addon_root/'inputs'/'4z80'/Path(triple_case['pdb']).name
    pdb.parent.mkdir()
    shutil.copy2(triple_case['pdb'],pdb)
    triple_case['pdb'] = str(pdb)
    shutil.copytree(root/'rb4'/'baseline'/'4z80',addon_root/'baseline'/'4z80')
    kernel.dump(addon_root/'plan.json',dict(rb=4,seed=42,arm='budget-forward',build=str(kernel.BUILD),
        baseline=str(BASE/'rb4'),cases=[triple_case],total_sequences=39,
        source_sha256={p.name:kernel.sha(p) for p in (addon_root/'source').iterdir()}))
    (addon_root/'READY').write_text('Authorized 4z80 RB4 optional-triple replacement; same frozen build and inputs.\n')
    for rb, cases in groups.items():
        group = root / f'rb{rb}'
        kernel.dump(group/'plan.json', dict(rb=rb, seed=42, arm='pair-only', build=str(kernel.BUILD),
            baseline=str(BASE/f'rb{rb}'), cases=cases, total_sequences=1103,
            sequence_validation='Reuse frozen RB1 current-build enumeration; identical source/cohort/PDB/sequence sets verified',
            enumeration_reference=str(RB1/'preflight'),
            source_sha256={p.name:kernel.sha(p) for p in (group/'source').iterdir()}))
        (group/'READY').write_text('Frozen RB1 build/runner, same-RB baselines, PDB hashes and sequence sets validated.\n')
    kernel.write_table(root/'resources.tsv', [{k:c[k] for k in ('rb','design','tier','cpus','gpus','gpu_type','mem_gib','heap_gib','host_gib','partition','time_limit','snapshot_availability','snapshot_candidate')} for c in all_cases])
    kernel.dump(root/'plan.json', dict(groups={str(rb):str(root/f'rb{rb}') for rb in groups},
        cases=all_cases, triple_addon=triple_case, total_sequences=2245, rb1_reference=str(RB1), build=str(kernel.BUILD)))
    # Validate all configurations before allowing any expensive task to start.
    from contextlib import ExitStack
    jobs = {2:[],4:[]}
    all_submissions = [dict(triple_case,arm='budget-forward'), *[dict(c,arm='pair-only') for c in all_cases]]
    with ExitStack() as stack:
        registries = {rb:stack.enter_context((root/f'rb{rb}'/'jobs.jsonl').open('x')) for rb in groups}
        registry = stack.enter_context((root/'jobs.jsonl').open('x'))
        addon_registry = stack.enter_context((addon_root/'jobs.jsonl').open('x'))
        for case in all_submissions:
            is_triple = case['arm']=='budget-forward'
            command = ['sbatch','--parsable','--account=grisman',
                '--partition='+case['partition'], '--job-name='+('triple4_4z80' if is_triple else 'pair38r'+str(case['rb'])+'_'+case['design']),
                '--cpus-per-task='+str(case['cpus']), '--mem='+str(case['mem_gib'])+'G',
                '--gres=gpu:'+case['gpu_type']+':'+str(case['gpus']), '--time='+case['time_limit'],
                '--export=ALL,PAIR38_ROOT='+case['group_root']+',CASE_INDEX='+str(case['index']),
                str(source/'run_pair38_rb24_20260922.slurm')]
            ident = subprocess.check_output(command, text=True).strip().split(';')[0]
            assert ident.isdigit(), ident
            record = dict(job=ident, rb=case['rb'], arm=case['arm'], index=case['index'], design=case['design'], command=command)
            for stream in (registry,addon_registry if is_triple else registries[case['rb']]):
                stream.write(json.dumps(record)+'\n'); stream.flush()
            if is_triple:
                triple_job = ident
                kernel.dump(addon_root/'submission.json',dict(jobs=[ident],systems=1,total_sequences=39))
            else:
                jobs[case['rb']].append(ident)
            print('SUBMITTED', ident, case['rb'], case['design'], case['tier'], flush=True)
    summaries = {}
    for rb, identifiers in jobs.items():
        summary = subprocess.check_output(['sbatch','--parsable','--account=grisman',
            '--job-name=pair38r'+str(rb)+'_summary', '--dependency=afterany:'+':'.join(identifiers),
            '--export=ALL,PAIR38_ROOT='+str(root/f'rb{rb}')+',PAIR38_MODE=summarize',
            str(source/'summarize_pair38_rb24_20260922.slurm')], text=True).strip().split(';')[0]
        summaries[rb] = summary
        kernel.dump(root/f'rb{rb}'/'submission.json',dict(jobs=identifiers,summary_job=summary,systems=38,total_sequences=1103))
    cross = subprocess.check_output(['sbatch','--parsable','--account=grisman',
        '--job-name=pair38_rb124_compare',
        '--dependency=afterany:'+':'.join([*summaries.values(),triple_job,'12680425']),
        '--export=ALL,PAIR38_ROOT='+str(root)+',PAIR38_MODE=compare',
        str(source/'summarize_pair38_rb24_20260922.slurm')], text=True).strip().split(';')[0]
    kernel.dump(root/'submission.json',dict(jobs=jobs, summaries=summaries, cross_rb_summary=cross,
        triple_addon_job=triple_job,total_compute_jobs=77,total_planned_sequences=2245,
        systems_per_rb=38,sequences_per_rb=1103, cpus=sum(c['cpus'] for c in all_submissions),
        gpus=sum(c['gpus'] for c in all_submissions), memory_gib=sum(c['mem_gib'] for c in all_submissions),
        tier_counts=dict(Counter(c['tier'] for c in all_submissions))))
    (root/'SUBMITTED').write_text('Both 38-system cohorts, 4z80 RB4 triple replacement, per-RB comparisons and cross-RB report submitted.\n')
    print('ROOT',root,'SUMMARY',summaries,'CROSS',cross,flush=True)


def compare(root):
    report = root/'cross_rb_comparison'
    report.mkdir(exist_ok=False)
    cohort = json.loads((RB1/'plan.json').read_text())['cases']
    totals, sequence_rows, graph_rows = [], [], []
    all_outcomes = {}
    for rb, group in ((1,RB1),(2,root/'rb2'),(4,root/'rb4')):
        planned = json.loads((group/'plan.json').read_text())
        pairs, bases = {}, {}
        completed = 0
        for case in planned['cases']:
            task_path = group/'tasks'/(str(case['index'])+'.json')
            task = json.loads(task_path.read_text()) if task_path.exists() else {}
            okay = task.get('technical_ok',False)
            # The independent final RB1 audit validated this known legacy
            # Unstable status case; do not edit its frozen task record.
            if rb==1 and case['design']=='3bu8' and task.get('status')=='AUDIT_FAILED':
                okay = task.get('java_exit')==0 and "'prot_status': 'Unstable'" in task.get('error','')
            completed += bool(okay)
            dest = Path(task['run']) if task.get('run') else None
            for arm, directory, outputs in [('pair',dest,pairs),('baseline',group/'baseline'/case['design'],bases)]:
                rr = kernel.rows(directory/(case['design']+'_packstar.csv')) if directory else []
                for r in rr:
                    if arm=='pair' and not okay:
                        continue
                    outputs[(case['design'],r['sequence'].strip())] = all(r[s+'_status']=='Estimated' for s in ('prot','lig','comp'))
            if dest and (dest/'run.log').exists():
                pattern = re.compile(r'InteractionGraph residual-budget cutoff: kept (\d+)/(\d+) edges \(cut (\d+)\), residualUpperBound=([\d.Ee+-]+), residualBudget=([\d.Ee+-]+)')
                event = 0
                with (dest/'run.log').open(errors='replace') as stream:
                    for line in stream:
                        m = pattern.search(line)
                        if m:
                            event += 1
                            kept, full, cut = map(int,m.group(1,2,3))
                            graph_rows.append(dict(rb=rb,design=case['design'],event=event,kept_edges=kept,
                                full_edges=full,cut_edges=cut,kept_fraction=kept/full if full else None,
                                residual_upper_bound=float(m[4]),residual_budget=float(m[5])))
        matched = set(pairs)&set(bases)
        totals.append(dict(rb=rb,complete_systems=completed,planned=1103,
            available_pair=len(pairs),pair_success=sum(pairs.values()),matched=len(matched),
            matched_pair_success=sum(pairs[k] for k in matched), matched_baseline_success=sum(bases[k] for k in matched),
            baseline_only_success=sum(bases[k] and not pairs[k] for k in matched),
            pair_only_success=sum(pairs[k] and not bases[k] for k in matched)))
        all_outcomes[rb] = dict(pair=pairs,baseline=bases)
    common = set.intersection(*(set(all_outcomes[rb]['pair'])&set(all_outcomes[rb]['baseline']) for rb in (1,2,4)))
    for row in totals:
        data = all_outcomes[row['rb']]
        row.update(common_all_rb=len(common),common_pair_success=sum(data['pair'][k] for k in common),
            common_baseline_success=sum(data['baseline'][k] for k in common))
    for case in cohort:
        for sequence in case['expected_sequences']:
            key=(case['design'],sequence)
            row=dict(design=key[0],sequence=key[1],in_common_all_rb=key in common)
            for rb in (1,2,4):
                for arm in ('pair','baseline'):
                    value=all_outcomes[rb][arm].get(key)
                    row[f'rb{rb}_{arm}']='Missing' if value is None else ('Estimated' if value else 'NonSuccess')
            sequence_rows.append(row)
    kernel.write_table(report/'summary.tsv',totals)
    kernel.write_table(report/'sequences.tsv',sequence_rows)
    if graph_rows:
        kernel.write_table(report/'graph_events.tsv',graph_rows)
    kernel.dump(report/'summary.json',totals)
    addon = root/'rb4_triple_4z80'
    addon_task_file = addon/'tasks/30.json'
    addon_task = json.loads(addon_task_file.read_text()) if addon_task_file.exists() else {}
    kernel.dump(report/'rb4_4z80_triple_task.json',addon_task)
    addon_rows = []
    if addon_task.get('technical_ok'):
        triple_rows = kernel.rows(Path(addon_task['run'])/'4z80_packstar.csv')
        for row in triple_rows:
            key=('4z80',row['sequence'].strip())
            pair=all_outcomes[4]['pair'].get(key)
            addon_rows.append(dict(design='4z80',sequence=key[1],
                pair_only='Missing' if pair is None else ('Estimated' if pair else 'NonSuccess'),
                triple_rerun='Estimated' if all(row[s+'_status']=='Estimated' for s in ('prot','lig','comp')) else 'NonSuccess',
                **{s+'_triple_status':row[s+'_status'] for s in ('prot','lig','comp')}))
        kernel.write_table(report/'rb4_4z80_pair_vs_triple_rerun.tsv',addon_rows)
    (report/'README.md').write_text('# Pair-only RB sensitivity\n\n'
        'summary.tsv reports each same-RB pair/triple comparison and the common sequence subset observed in all six arms. '
        'Success requires all three partition functions Estimated. Missing/process failures are separate from NonSuccess. '
        'sequences.tsv retains all 1103 planned sequence identities. graph_events.tsv records actual retained/cut edges; '
        'events may repeat cached partition functions and are not independent samples. '
        'RB1 3bua baseline is missing after a host-budget error; RB4 4z80 baseline is missing after a two-day timeout. '
        'The separately authorized 4z80 RB4 optional-triple replacement is compared in rb4_4z80_pair_vs_triple_rerun.tsv, '
        'with its technical status in rb4_4z80_triple_task.json; it does not overwrite the historical baseline. '
        'This is a seed42 historical comparison with differing hardware/root budgets; timing is not a speedup metric.\n')
    print('CROSS_RB_REPORT',report,flush=True)


if __name__=='__main__':
    assert os.environ.get('SLURM_JOB_ID')
    assert os.environ.get('SLURM_JOB_ACCOUNT','grisman')=='grisman'
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['launch','run','summarize','compare'])
    parser.add_argument('--root',type=Path)
    parser.add_argument('--index',type=int)
    args=parser.parse_args()
    if args.root:
        assert args.root.absolute().is_relative_to(PARENT)
    if args.mode=='launch': launch()
    elif args.mode=='run': kernel.run(args.root,args.index)
    elif args.mode=='summarize': kernel.summarize(args.root)
    else: compare(args.root)
