"""Export the requested small versioned plan snapshot from a Slurm report."""
import csv
import json
import os
from pathlib import Path

assert os.environ.get('SLURM_JOB_ID'), 'Run the report export through Slurm.'
here=Path(__file__).resolve().parent
report=(here/'latest_report').resolve()
summary=json.loads((report/'table_summary.json').read_text())
supplement=json.loads((here.parent/'plan_supplements_20261007/registry.json').read_text())
preferences=json.loads((here/'reporting_preferences.json').read_text())
with (report/'unified_plan.tsv').open(newline='') as handle:
    rows=list(csv.DictReader(handle,delimiter='\t'))
fields=['system','planned_tier_count','design','state','job','destination','MARK_display','current_MARK_CPU',
        'PACK_display','new_PACK_job','history_display','recommendation','dependency']
with (here/'current_plan.tsv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
    writer.writeheader()
    writer.writerows(dict(r,dependency=r.get('dependency') or 'NA') for r in rows)
text='# MARK* 统一计划表\n\n'
text+=f'快照：{summary["snapshot"]}。共 {summary["systems"]} 个系统、{summary["planned_tiers"]} 档，其中 {summary["completed_results_retained"]} 档已有本轮 MARK* 结果。\n\n'
text+='计划档数包含条件保留档。只有本轮 MARK* 已出结果的单元格使用斜体。G 表示四张 PRO 6000 GPU，C 表示 CPU 批次。★ 为历史完成档中最接近14天的一档。完成表示有完整结果表，不保证每条序列均达到 Estimated。\n\n'
text+='“暂不排队”的三个条件保留档已撤销未启动提交；“待提交”仍是计划项。“DCC计划”尚无远端提交回执，“DCC待核”表示曾提交、当前进度未核实。\n\n'
limits=preferences.get('concurrency_limits')
if limits:
    text+=f'本地并发上限：compsci {limits["compsci"]} + fennario {limits["fennario"]} = {limits["compsci"]+limits["fennario"]}，每任务 64 CPU、192 GiB 内存；实际运行数以本表快照为准。DCC 计划上限仍为 {limits["dcc"]}，远端状态未核实。详见[并发调整记录](../local_concurrency14_20261007/registry.json)。\n\n'
text+=f'PACK* 补测任务 `{supplement["pack_job"]}` 包含 `5dc0_p9`、`5d68_p4`、`3k3q_m1`，使用四张 PRO 6000 GPU。\n\n'
text+='[逐档任务号与依赖](current_plan.tsv) · [本地队列调整记录](../queue_alignment_20261007/registry.json) · [DCC 提交清单](../deadline_focus_20261006/dcc_launch_priority.tsv)\n\n'
text+=(report/'unified_table.txt').read_text()
text+='\n## 5dc0 PACK* 记录\n\n'+(report/'5dc0_table.txt').read_text()
(here/'CURRENT_PLAN.md').write_text(text)
assert len(rows)==summary['planned_tiers']
print('EXPORTED',here/'CURRENT_PLAN.md',flush=True)
