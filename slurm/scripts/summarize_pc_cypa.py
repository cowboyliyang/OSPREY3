#!/usr/bin/env python3
"""Summarize a completed CypA run without treating NMR mapping as proven."""
import csv
import json
import math
import os
from pathlib import Path
import statistics
import sys

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Slurm required')
root = Path(sys.argv[1])
out = root / 'output'
if not (out / 'COMPUTATION_FINISHED').exists():
    raise SystemExit('Computation has not finished; no success report generated')
with (out / 'results.tsv').open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
if len(rows) != 12:
    raise SystemExit('Expected all four variants and all three seeds')
summary = []
for variant in ('WT','S99T','S99T_C115S','S99T_C115S_I97V'):
    group = [r for r in rows if r['variant'] == variant]
    if {r['seed'] for r in group} != {'42','43','44'}:
        raise SystemExit('Missing or duplicated seed')
    valid = [float(r['P_C']) for r in group if r['pc_status'] == 'RESOLVED']
    if any(not math.isfinite(p) or p < 0 or p > 1 for p in valid):
        raise SystemExit('Invalid resolved probability')
    ref = float(group[0]['P_ref'])
    nmr = float(group[0]['NMR_pB'])
    item = dict(variant=variant, P_ref=ref, resolved_seeds=len(valid),
                PC_min=min(valid) if valid else None, PC_max=max(valid) if valid else None,
                PC_mean=statistics.mean(valid) if valid else None,
                NMR_pB=nmr if math.isfinite(nmr) else None,
                NMR_reported_error=float(group[0]['NMR_reported_error']) if math.isfinite(nmr) else None,
                within_model_all_pass=all(r['within_model_check']=='PASS_ABS_0.02' for r in group),
                P_ref_abs_difference_NMR=abs(ref-nmr) if math.isfinite(nmr) else None,
                event_sensitivity_min=float(group[0]['P_ref_30_90']),
                event_sensitivity_max=float(group[0]['P_ref_m15_135']))
    summary.append(item)
report = dict(computation_finished=True, seed_rows=len(rows), variants=summary,
              within_model_all_pass=all(s['within_model_all_pass'] for s in summary),
              biological_ground_truth_validation='NOT_ESTABLISHED',
              reason='Local Phe113-in predicate and finite fixed-backbone measure are approximations to NMR group-I B state; no event confidence intervals.',
              preparation_qc='separate structural audit required')
(out / 'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
lines = ['# CypA external comparison','',
         '| Variant | P_ref | PACK* P_C mean (resolved seeds) | NMR pB | Model check |',
         '|---|---:|---:|---:|---|']
for r in summary:
    pc = 'unresolved' if r['PC_mean'] is None else '{:.6g} ({}/3)'.format(r['PC_mean'],r['resolved_seeds'])
    nmr = 'not specified' if r['NMR_pB'] is None else str(r['NMR_pB'])
    lines.append('| {} | {:.6g} | {} | {} | {} |'.format(r['variant'],r['P_ref'],pc,nmr,'pass' if r['within_model_all_pass'] else 'not passed'))
lines += ['',report['reason'], '', 'This is an executed external comparison; completion does not establish biological accuracy.']
(out / 'summary.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(report,indent=2,allow_nan=False))
