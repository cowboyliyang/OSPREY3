"""Audit latest pair-only final batches. Must execute through Slurm."""
import csv
import datetime
import json
import math
import os
from pathlib import Path
import re
import subprocess

ROOT = Path('/usr/xtmp/lz280/packstar_pair38_rb24_20260922/launch_12681148')
RT = 0.0019891 * 298.15


def read_rows(path, delimiter='\t'):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def keys(path):
    return {r['key']: r['value'] for r in read_rows(path)} if path.exists() else {}


def write_table(path, data):
    if not data:
        path.write_text('')
        return
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(data)


def verify_samples(path, n, k, severity):
    rows = read_rows(path)
    assert len(rows) == n, path
    logs = []
    for row in rows:
        relative = float(row['logRelativeWeight']) - 1 / RT
        assert (relative > 0) == (row['tailEvent'] == 'true'), path
        if relative > 0:
            logs.append(relative + math.log1p(-math.exp(-relative)) if relative > 40 else math.log(math.expm1(relative)))
    assert len(logs) == k and k > 0, path
    high = max(logs)
    log_mean = high + math.log(math.fsum(math.exp(x-high) for x in logs)) - math.log(k)
    if math.isfinite(severity):
        assert math.isclose(log_mean, math.log(severity), abs_tol=1e-9, rel_tol=1e-10), path
    else:
        assert log_mean > 709 or high > 709, path
    return log_mean / math.log(10), len({r['conf'] for r in rows if r['tailEvent'] == 'true'})


def main():
    job = os.environ['SLURM_JOB_ID']
    assert os.environ.get('SLURM_JOB_ACCOUNT', 'grisman') == 'grisman'
    out = Path('/usr/xtmp/lz280') / f'packstar_rb24_final_s0_audit_20260924_{job}'
    out.mkdir(exist_ok=False)
    records, systems, anomalies = [], [], []
    for rb in (2, 4):
        group = ROOT / f'rb{rb}'
        plan = json.loads((group/'plan.json').read_text())
        assert plan['rb'] == rb and plan['arm'] == 'pair-only' and plan['seed'] == 42
        assert len(plan['cases']) == 38
        for case in sorted(plan['cases'], key=lambda c:c['design']):
            design = case['design']
            runs = list((group/'runs').glob(design+'_full_pair-only_s42_*'))
            assert len(runs) == 1, (design, runs)
            run = runs[0]
            m = json.loads((run/'manifest.json').read_text())
            p = m['properties']
            assert p['branchdp.cutoff.residualBudget'] == str(rb)
            assert p['packstar.pac.frequencySeverity.tripleEta'] == 'false'
            assert float(p['packstar.pac.frequencySeverity.severityCap']) == 20
            assert float(p['packstar.pac.frequencySeverity.relativeBoundKcal']) == 1
            result_path = run/(design+'_packstar.csv')
            result_rows = read_rows(result_path, ',') if result_path.exists() else []
            # The log identifies the physical state without matching rounded endpoints.
            roles = {}
            role = None
            for line in (run/'run.log').read_text(errors='replace').splitlines():
                match = re.search(r'estimator activated for state=(\w+)', line)
                if match:
                    role = match[1].lower()
                match = re.search(r'identityHash=([0-9a-f]+)', line)
                if match:
                    roles['state-'+match[1]] = role
            states = sorted((run/'adaptive_frequency_severity').glob('state-*'))
            local = []
            partial = 0
            for state in states:
                path = state/'frequency_severity_final.tsv'
                if not path.exists():
                    if (state/'frequency_severity_final_samples.tsv').exists():
                        partial += 1
                    continue
                rows = read_rows(path)
                assert len(rows) == 1 and rows[0]['stage'] == 'final', path
                row = rows[0]
                n, k = int(row['n']), int(row['tailCount'])
                s0 = float(row['empiricalConditionalSeverity'])
                assert 0 <= k <= n and not math.isnan(s0) and s0 >= 0
                if k:
                    if math.isfinite(s0):
                        assert math.isclose(s0, float(row['empiricalTailMean'])*n/k, rel_tol=1e-9, abs_tol=1e-9)
                else:
                    assert s0 == 0
                final = keys(state/'frequency_severity_final_interval.tsv')
                fail = keys(state/'frequency_severity_failure.tsv')
                record = dict(rb=rb, design=design, state_id=state.name, role=roles.get(state.name,''),
                    n=n, tail_count=k, observed_s0=s0 if k else None,
                    max_single_excess=float(row['observedMaxConditionalSeverity']),
                    over_20=bool(k and s0>20), certificate_valid=final.get('certificateValid',''),
                    target_reached=final.get('targetReached',''), final_complete=bool(final),
                    terminal_failure=fail.get('reason',''), run_status=m['status'],
                    sequence='', result_status='', verified_log10_s0='', tail_unique_confs='', path=str(path))
                if record['over_20']:
                    record['verified_log10_s0'], record['tail_unique_confs'] = verify_samples(
                        state/'frequency_severity_final_samples.tsv', n, k, s0)
                    matches = []
                    if final and 'logZLower' in final and 'logZUpper' in final:
                        lo, hi = (float(final[key])/math.log(10) for key in ('logZLower','logZUpper'))
                        for result in result_rows:
                            for short, name in [('prot','protein'),('lig','ligand'),('comp','complex')]:
                                if record['role'] and name != record['role']:
                                    continue
                                if (abs(float(result[short+'_qstar_lb_log10'])-lo)<1e-5 and
                                    abs(float(result[short+'_qstar_ub_log10'])-hi)<1e-5):
                                    matches.append((name,result['sequence'],result[short+'_status']))
                    if matches:
                        assert len({x[0] for x in matches}) == 1
                        record['role'] = matches[0][0]
                        record['sequence'] = '; '.join(sorted({x[1] for x in matches}))
                        record['result_status'] = '; '.join(sorted({x[2] for x in matches}))
                    else:
                        anomalies.append(dict(rb=rb, design=design, state=state.name, issue='no CSV interval match'))
                records.append(record)
                local.append(record)
            systems.append(dict(rb=rb,design=design,status=m['status'],exit_code=m.get('exit_code'),
                ended='end_epoch' in m,result_rows=len(result_rows),expected_rows=case['expected_rows'],
                states=len(states),final_batches=len(local),partial_final_artifacts=partial,
                final_with_tail=sum(r['tail_count']>0 for r in local),
                final_over_20=sum(r['over_20'] for r in local),run=str(run)))
    summary = dict(root=str(ROOT),job=job,as_of=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        definition='Per distinct final state batch: mean(R/C - 1 among R>C) > 20. No-tail batches have undefined severity.',
        deduplication='One case per run/state directory; cached ligand/protein reused across sequence CSV rows counted once.',
        groups={},anomalies=anomalies,output=str(out))
    for rb in (2,4):
        rr = [r for r in records if r['rb']==rb]
        ss = [s for s in systems if s['rb']==rb]
        xx = [r for r in rr if r['over_20']]
        summary['groups'][str(rb)] = dict(systems=len(ss),systems_ended=sum(s['ended'] for s in ss),
            unfinished_systems=[s['design'] for s in ss if not s['ended']],
            final_batches=len(rr),with_tail=sum(r['tail_count']>0 for r in rr),
            over_20=len(xx),over_20_estimated=sum(r['result_status']=='Estimated' for r in xx),
            systems_over_20=sorted({r['design'] for r in xx}),
            max_single_over_20=sum(r['max_single_excess']>20 for r in rr),
            partial_final_artifacts=sum(s['partial_final_artifacts'] for s in ss))
    exceeds = sorted((r for r in records if r['over_20']), key=lambda r:(r['rb'],-r['observed_s0']))
    write_table(out/'all_final.tsv',records)
    write_table(out/'exceeds_20.tsv',exceeds)
    write_table(out/'per_system.tsv',systems)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Latest pair-only rb=2,4: final conditional severity audit','',summary['definition'],'',
        '| rb | Final batches | With tail | Mean severity >20 | Estimated among >20 | Systems with >20 |',
        '|---|---:|---:|---:|---:|---|']
    for rb,g in summary['groups'].items():
        lines.append(f"| {rb} | {g['final_batches']} | {g['with_tail']} | {g['over_20']} | {g['over_20_estimated']} | {', '.join(g['systems_over_20'])} |")
    lines += ['', '| rb | System | Sequence | State | Observed severity | Tail / N | Result status |',
        '|---|---|---|---|---:|---:|---|']
    for r in exceeds:
        lines.append(f"| {r['rb']} | {r['design']} | {r['sequence'] or r['state_id']} | {r['role']} | {r['observed_s0']:.9g} | {r['tail_count']}/{r['n']} | {r['result_status'] or 'unmapped'} |")
    lines += ['', 'As of '+summary['as_of'], '', 'Unfinished designs: '+json.dumps({rb:g['unfinished_systems'] for rb,g in summary['groups'].items()}),
        '', 'Raw final samples independently recomputed for every case above 20.', '', 'Anomalies: '+json.dumps(anomalies)]
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    print((out/'report.md').read_text(),flush=True)
    print('OUTPUT='+str(out),flush=True)


if __name__ == '__main__':
    main()
