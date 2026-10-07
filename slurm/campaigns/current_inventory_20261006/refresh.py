"""Collect current MARK* completion artifacts, live local jobs and DCC plans."""
import csv
import json
import os
import re
import subprocess
from collections import Counter
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

assert os.environ.get('SLURM_JOB_ID'), 'Submit this inventory refresh through Slurm.'
campaign = Path(__file__).resolve().parent
campaigns = campaign.parent
old = Path('/usr/xtmp/lz280/current_markstar_tiers_12827815')
out = Path('/usr/xtmp/lz280/markstar_current_inventory_' + os.environ['SLURM_JOB_ID'])
out.mkdir(parents=True, exist_ok=False)

def tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def duration_seconds(text):
    days = 0
    if '-' in text:
        days, text = text.split('-', 1)
        days = int(days)
    parts = [int(x) for x in text.split(':')]
    while len(parts) < 3:
        parts.insert(0, 0)
    return days * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]

def estimated(row):
    return all(row.get(part + '_status') == 'Estimated' for part in ['prot', 'lig', 'comp'])

def clock_text(seconds):
    seconds = round(seconds)
    return f'{seconds//3600:02d}:{seconds%3600//60:02d}:{seconds%60:02d}'

exclusions_path = campaign / 'excluded_designs.json'
excluded = set(json.loads(exclusions_path.read_text())['designs']) if exclusions_path.exists() else set()
definitions = {r['design']: r for r in tsv(campaigns / 'deadline_focus_20261006/design_actions.tsv') if r['design'] not in excluded}
rows = tsv(old / 'all_tiers.tsv')
columns = list(rows[0])
held = {key: '' for key in columns}
held.update(design='2xxm_flex_p3', system='2xxm', tier='3', state='备用暂停', job='12812439_43',
            destination='compsci', historical_MARK='38.22h', PACK_minutes_hardware='本轮测量待核')
rows.append(held)
rows = [r for r in rows if r['design'] not in excluded]
queue_raw = subprocess.check_output(['squeue', '-h', '-r', '-u', 'lz280', '-o', '%i|%j|%T|%M|%R|%P|%E'], text=True)
(out / 'queue.txt').write_text(queue_raw)
queue = {}
for line in queue_raw.splitlines():
    fields = line.split('|', 6)
    if len(fields) == 7:
        queue[fields[0]] = dict(zip(['job', 'name', 'state', 'elapsed', 'reason', 'partition', 'dependency'], fields))
dcc_priority = {r['design_id']: r for r in tsv(campaigns / 'deadline_focus_20261006/dcc_launch_priority.tsv')}
known = {r['design'] for r in rows}
for definition in tsv(campaigns / 'deadline_focus_20261006/dcc_new_designs.tsv'):
    name = definition['design_id']
    if name not in known and name in dcc_priority and name not in excluded:
        row = {key: '' for key in columns}
        row.update(design=name, system=definition['system'], tier=definition['flex_delta'],
                   state='DCC新增计划，无提交回执', destination='dcc', historical_MARK='无')
        rows.append(row)
        known.add(name)
conditional_path = campaign/'conditional_deferred.json'
conditional = json.loads(conditional_path.read_text())['designs'] if conditional_path.exists() else {}
dcc_complete = {r['design_id']: r for r in tsv(Path('/usr/xtmp/lz280/markstar_dcc64_completed38_20261006/design_summary.tsv'))}
local_root = Path('/usr/xtmp/lz280/markstar_local57_cpu64_20261005/A12812439')
for row in rows:
    name = row['design']
    row.update(MARK_hms='', MARK_elapsed_hms='', dependency='', MARK_fully_estimated='',
               MARK_rows='', MARK_expected_rows='', MARK_output='', provenance='', PACK_status='',
               historical_local_copy_job='')
    if name in dcc_complete:
        record = dcc_complete[name]
        # TRANSFER uses its own header labels; use the archived CSV for row quality.
        directory = Path('/usr/xtmp/lz280/markstar_dcc64_completed38_20261006/designs') / name
        result_rows = list(csv.DictReader((directory / (name + '_markstar.csv')).open()))
        seconds = float(record['elapsed_seconds'])
        row.update(state='完成', MARK_output=str(directory), MARK_rows=len(result_rows),
                   MARK_fully_estimated=sum(estimated(r) for r in result_rows),
                   MARK_expected_rows=int(record['expected_rows']), MARK_hms=clock_text(seconds),
                   current_MARK_hours=round(seconds/3600, 6), current_MARK_CPU=record['cpu_model'],
                   job=record['array_job_id'] + '_' + record['array_task_id'],
                   provenance='DCC completed38 export; remote live state unavailable')
        if row['PACK_seconds']:
            row['raw_wall_ratio'] = round(seconds/float(row['PACK_seconds']), 4)
    root = local_root / name
    manifest_path = root / 'run_manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        row['MARK_output'] = str(root)
        row['provenance'] = 'Local run artifact plus current scheduler query'
        csv_path = root / (name + '_markstar.csv')
        if manifest['status'] == 'COMPLETED' and csv_path.exists():
            result_rows = list(csv.DictReader(csv_path.open()))
            seconds = manifest['end_epoch'] - manifest['start_epoch']
            wall_path = root / 'wall.time'
            if wall_path.exists():
                found = re.search(r'elapsed=([0-9.]+)', wall_path.read_text())
                if found:
                    seconds = float(found[1])
            row.update(state='完成', current_MARK_hours=round(seconds/3600, 6),
                       MARK_hms=clock_text(seconds), MARK_rows=len(result_rows),
                       MARK_expected_rows=len(result_rows), MARK_fully_estimated=sum(estimated(r) for r in result_rows))
            cpu = root / 'lscpu.txt'
            if cpu.exists():
                found = re.search(r'^Model name:\s*(.+)$', cpu.read_text(), flags=re.M)
                if found:
                    row['current_MARK_CPU'] = found[1].strip()
            if row['PACK_seconds']:
                row['raw_wall_ratio'] = round(seconds/float(row['PACK_seconds']), 4)
                pack_root = Path(row['PACK_output'])
                pack_csv = pack_root / (name + '_packstar.csv')
                if pack_csv.exists():
                    pack_rows = {r['sequence']: r for r in csv.DictReader(pack_csv.open())}
                    matched = sum(estimated(r) and r['sequence'] in pack_rows and estimated(pack_rows[r['sequence']]) for r in result_rows)
                    row['joint_verified'] = f'{matched}/{len(result_rows)}'
    if row['destination'] == 'dcc' and row['state'] != '完成':
        if row['job']:
            row['historical_local_copy_job'] = row['job']
            row['job'] = ''
        prior = dcc_priority[name]['prior_submission']
        row['state'] = 'DCC原提交，进度待核' if prior.startswith('ALREADY_SUBMITTED') else 'DCC新增计划，无提交回执'
        row['provenance'] = 'DCC handoff manifest; no current remote scheduler snapshot'
    elif row['job'] in queue:
        live = queue[row['job']]
        row['state'] = '运行' if live['state'] == 'RUNNING' else '备用暂停' if live['reason'] == '(JobHeldUser)' else '排队'
        row['MARK_elapsed_hms'] = clock_text(duration_seconds(live['elapsed']))
        row['dependency'] = live['dependency']
        row['destination'] = 'compsci' if live['partition'] == 'compsci' else 'fennario'
    if name in conditional and row['state'] not in ('完成','运行') and row['job'] not in queue:
        row.update(state='条件保留，暂不排队', job='', dependency='', MARK_elapsed_hms='',
                   provenance='User deferred this tier; prior pending task cancelled; '+conditional[name]['audit'])

# Existing speed ratios are raw whole-run ratios; the matched Estimated count is separate.
by_name = {r['design']: r for r in rows}
pack_groups = ['packstar_gpu101_20261005/run_12814879', 'markstar_near14_20261006/pack_12827800', 'markstar_deadline_focus_20261006/pack_12827802']
supplement_path = campaigns/'plan_supplements_20261007/registry.json'
if supplement_path.exists():
    supplement = json.loads(supplement_path.read_text())
    if supplement.get('pack_output'):
        pack_groups.append(supplement['pack_output'])
for group in pack_groups:
    status_path = Path('/usr/xtmp/lz280') / group / 'status.tsv'
    if not status_path.exists():
        continue
    for record in tsv(status_path):
        if record['design'] not in by_name:
            continue
        row = by_name[record['design']]
        row['PACK_status'] = record['status']
        if record['status'] not in ['COMPLETED', 'INCOMPLETE_ESTIMATES'] or not record['output']:
            continue
        path = Path(record['output'])
        wall = path / 'wall.time'
        if wall.exists():
            found = re.search(r'elapsed=([0-9.]+)', wall.read_text())
            if found:
                seconds = float(found[1])
                row['PACK_seconds'] = seconds
                row['PACK_minutes_hardware'] = f'{seconds/60:.2f}G'
                row['PACK_output'] = str(path)
for row in rows:
    if row['state'] == '完成' and not row['MARK_hms'] and row['current_MARK_hours']:
        row['MARK_hms'] = clock_text(float(row['current_MARK_hours']) * 3600)

systems = sorted({r['system'] for r in rows})
by_system = {s: sorted([r for r in rows if r['system'] == s], key=lambda r: int(r['tier'])) for s in systems}
system_rows = []
for system in systems:
    members = by_system[system]
    def tiers(state):
        return '、'.join(r['design'].split('_flex_')[1] for r in members if r['state'] == state) or '—'
    completed = sum(r['state'] == '完成' for r in members)
    system_rows.append(dict(system=system, completed=tiers('完成'), running=tiers('运行'),
                            local_pending=tiers('排队'), dcc_previous_submission_unverified=tiers('DCC原提交，进度待核'),
                            dcc_new_plan_without_receipt=tiers('DCC新增计划，无提交回执'), held=tiers('备用暂停'),
                            completed_tier_count=completed, more_completed_tiers_needed_by_oct16=max(0, 2-completed)))

states = Counter(r['state'] for r in rows)
snapshot = datetime.now(ZoneInfo('America/New_York')).isoformat()
summary = dict(snapshot=snapshot, systems=len(systems), designs_in_current_inventory=len(rows), states=dict(states),
               completed_systems=sum(r['completed_tier_count'] > 0 for r in system_rows),
               systems_with_two_completed_tiers=sum(r['completed_tier_count'] >= 2 for r in system_rows),
               additional_completed_tiers_needed_for_two_per_system=sum(r['more_completed_tiers_needed_by_oct16'] for r in system_rows),
               goal='Two completed current-round tiers per system by October 16; final 14-day round starts October 16.',
               status_scope='Main 64-CPU MARK* campaign. Completed means normal run completion with result rows, not all rows Estimated.',
               speedup_scope='Raw whole-run wall ratios; joint_verified records jointly Estimated rows only and is not a numerical agreement test.',
               dcc_scope='38 returned completions verified from transfer. Original remote submissions and new handoff definitions are distinct; current remote progress is unavailable.',
               legacy_48cpu_comparisons=[queue[j] for j in ['12825151','12825152'] if j in queue],
               output=str(out), source_snapshot=str(old), coverage_rebalance_registry=str(campaigns / 'coverage_lane_rebalance_20261006/registry.json'))
write_tsv(out / 'designs.tsv', sorted(rows, key=lambda r: (r['system'], int(r['tier']))))
write_tsv(out / 'systems.tsv', system_rows)
active = {r['design']: r for r in rows}
definition_rows = []
for name, definition in sorted(definitions.items()):
    current = active.get(name)
    definition_rows.append(dict(design=name, system=definition['system'], policy_action=definition['action'],
                                current_state=current['state'] if current else '当前计划暂缓；DCC实际队列未核实' if definition['current_assignment'] == 'dcc' else '当前计划暂缓',
                                job=current['job'] if current else '', historical_MARK_hours=definition['old_mark_hours'],
                                historical_PACK_minutes=definition['old_pack_minutes']))
write_tsv(out / 'all_definitions.tsv', definition_rows)
summary['all_definition_count'] = len(definition_rows)
summary['excluded_designs'] = sorted(excluded)
summary['exclusions_manifest'] = str(exclusions_path)
(out / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')

lines = [f'# 当前 MARK* design 汇总\n\n快照：{snapshot}。\n',
         '完成表示任务正常结束并输出结果；并非所有序列均达到 Estimated。DCC 实时状态未核实。\n',
         '| 系统 | 已完成 | 正在跑 | 本地排队 | DCC 原提交，进度待核 | DCC 新增计划，无提交回执 | 备用暂停 |',
         '| --- | --- | --- | --- | --- | --- | --- |']
for r in system_rows:
    lines.append('| ' + ' | '.join(str(r[key]) for key in ['system','completed','running','local_pending','dcc_previous_submission_unverified','dcc_new_plan_without_receipt','held']) + ' |')
lines.extend(['\n全名为“系统_flex_档位”。\n',
              f'已完成 design：{states["完成"]}；本地运行：{states["运行"]}；本地排队：{states["排队"]}；备用暂停：{states["备用暂停"]}。\n',
              f'已有两档完成结果的系统：{summary["systems_with_two_completed_tiers"]}/38；还需增加 {summary["additional_completed_tiers_needed_for_two_per_system"]} 个完成档位。\n',
              '[逐 design 明细](designs.tsv)：含本轮 MARK* 耗时、PACK* 参考耗时、历史 MARK* 耗时、CPU、job ID、当前依赖和结果路径。\n',
              '[全部定义](all_definitions.tsv)：同时保留当前暂缓的档位。\n'])
(out / 'overview.md').write_text('\n'.join(lines)+'\n')
write_tsv(campaign / 'systems.tsv', system_rows)
(campaign / 'registry.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')
print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
