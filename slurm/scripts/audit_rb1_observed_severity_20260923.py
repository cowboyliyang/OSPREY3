"""Read-only audit of observed conditional severity; execute through Slurm."""
import csv
import json
import math
import os
from pathlib import Path
import re
import statistics

ROOT = Path('/usr/xtmp/lz280/packstar_pair38_20260921/launch_12677939')
RT = (1.9891 / 1000.0) * 298.15


def verify_raw(path, expected_n, expected_count, expected_severity):
    samples = rows(path.with_name(path.stem + '_samples.tsv'))
    assert len(samples) == expected_n
    log_excesses = []
    for sample in samples:
        relative = float(sample['logRelativeWeight']) - 1.0 / RT
        assert (relative > 0) == (sample['tailEvent'] == 'true')
        if relative > 0:
            log_excesses.append(relative if relative > 50 else math.log(math.expm1(relative)))
    assert len(log_excesses) == expected_count
    high = max(log_excesses)
    log_mean = high + math.log(math.fsum(math.exp(x-high) for x in log_excesses)) - math.log(expected_count)
    if math.isfinite(expected_severity):
        assert math.isclose(log_mean, math.log(expected_severity), rel_tol=1e-10, abs_tol=1e-9), path
    else:
        assert log_mean > 709 or high > 709, path
    return log_mean / math.log(10)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def keys(path):
    return {r['key']: r['value'] for r in rows(path)} if path.exists() else {}


def table(path, records):
    if not records:
        path.write_text('')
        return
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def describe(records):
    positive = [r for r in records if r['tail_count'] > 0]
    values = sorted(r['observed_s0'] for r in positive)
    return dict(batches=len(records), draws=sum(r['n'] for r in records),
                tail_observed_batches=len(positive),
                no_tail_batches=len(records)-len(positive),
                tail_draws=sum(r['tail_count'] for r in positive),
                tail_batch_median=statistics.median(values) if values else None,
                tail_batch_max=max(values) if values else None,
                batches_over_20=sum(v > 20 for v in values),
                systems_over_20=sorted({r['design'] for r in positive if r['observed_s0'] > 20}))


def main():
    job = os.environ['SLURM_JOB_ID']
    assert os.environ.get('SLURM_JOB_ACCOUNT', 'grisman') == 'grisman'
    out = Path('/usr/xtmp/lz280') / ('packstar_rb1_s0_audit_20260923_' + job)
    out.mkdir(exist_ok=False)
    plan = json.loads((ROOT/'plan.json').read_text())
    assert plan['rb'] == 1 and plan['arm'] == 'pair-only' and len(plan['cases']) == 38
    inventory, artifacts, systems, errors = [], [], [], []
    for case in sorted(plan['cases'], key=lambda c:c['design']):
        design = case['design']
        run_dirs = list((ROOT/'runs').glob(design + '_full_pair-only_s42_*'))
        assert len(run_dirs) == 1, (design, run_dirs)
        run = run_dirs[0]
        manifest = json.loads((run/'manifest.json').read_text())
        assert manifest['properties']['branchdp.cutoff.residualBudget'] == '1'
        assert manifest['properties']['packstar.pac.frequencySeverity.tripleEta'] == 'false'
        states = sorted((run/'adaptive_frequency_severity').glob('state-*'))
        state_status = {}
        for state in states:
            protocol = keys(state/'frequency_severity_protocol.tsv')
            final = keys(state/'frequency_severity_final_interval.tsv')
            failure = keys(state/'frequency_severity_failure.tsv')
            state_status[state.name] = dict(final=bool(final), valid=final.get('certificateValid'),
                                            failure=failure.get('reason', ''))
            if protocol:
                assert float(protocol['severityCapS0']) == 20
                assert float(protocol['relativeBoundKcal']) == 1
            for path in sorted(state.glob('frequency_severity_*.tsv')):
                if path.name.endswith('_samples.tsv'):
                    continue
                with path.open() as stream:
                    reader = csv.DictReader(stream, delimiter='\t')
                    if 'empiricalConditionalSeverity' not in (reader.fieldnames or []):
                        continue
                    batch = list(reader)
                assert len(batch) == 1, path
                row = batch[0]
                inventory.append(dict(path=str(path), bytes=path.stat().st_size))
                n, count = int(row['n']), int(row['tailCount'])
                severity = float(row['empiricalConditionalSeverity'])
                assert not math.isnan(severity) and severity >= 0, (path, severity)
                assert 0 <= count <= n and float(row['severityCapS0']) == 20
                if count and math.isfinite(severity):
                    assert math.isclose(severity, float(row['empiricalTailMean'])*n/count,
                                        rel_tol=1e-10, abs_tol=1e-10)
                elif not count:
                    assert severity == 0
                stage = row['stage']
                verified_log10 = verify_raw(path, n, count, severity) if severity > 20 else None
                matches = []
                if severity > 20 and final:
                    lo, hi = float(final['logZLower'])/math.log(10), float(final['logZUpper'])/math.log(10)
                    with (run/(design+'_packstar.csv')).open() as stream:
                        for result in csv.DictReader(stream):
                            for role in ('prot', 'lig', 'comp'):
                                if (abs(float(result[role+'_qstar_lb_log10'])-lo)<1e-5
                                        and abs(float(result[role+'_qstar_ub_log10'])-hi)<1e-5):
                                    matches.append(role+':'+result['sequence']+':'+result[role+'_status'])
                artifacts.append(dict(design=design, state=state.name, stage=stage,
                    phase='final' if stage == 'final' else 'prefinal',
                    n=n, tail_count=count, observed_s0=severity if count else None,
                    max_single_excess=float(row['observedMaxConditionalSeverity']),
                    severity_test_rejected=row['severityTestRejected'],
                    certificate_valid=final.get('certificateValid', ''),
                    terminal_failure=failure.get('reason', ''),
                    verified_log10_s0=verified_log10, final_result_matches=';'.join(matches), path=str(path)))
        local = [r for r in artifacts if r['design'] == design]
        final_rows = [r for r in local if r['phase'] == 'final']
        pre_rows = [r for r in local if r['phase'] == 'prefinal']
        all_stats, final_stats, pre_stats = map(describe, [local, final_rows, pre_rows])
        worst = max((r for r in local if r['tail_count']), key=lambda r:r['observed_s0'], default=None)
        systems.append(dict(design=design, states=len(states), final_batches=final_stats['batches'],
            final_tail_batches=final_stats['tail_observed_batches'], final_tail_draws=final_stats['tail_draws'],
            final_max_s0=final_stats['tail_batch_max'], prefinal_max_s0=pre_stats['tail_batch_max'],
            all_max_s0=all_stats['tail_batch_max'], batches_over_20=all_stats['batches_over_20'],
            final_batches_over_20=final_stats['batches_over_20'],
            all_tail_batches=all_stats['tail_observed_batches'],
            worst_stage=worst['stage'] if worst else '', worst_tail_count=worst['tail_count'] if worst else 0,
            worst_state=worst['state'] if worst else '', run=str(run)))
        print('AUDITED', design, json.dumps(systems[-1]), flush=True)
    exceeds = sorted((r for r in artifacts if r['tail_count'] and r['observed_s0'] > 20),
                     key=lambda r:r['observed_s0'], reverse=True)
    summary = dict(root=str(ROOT), job=job, systems=len(systems),
        definition='mean(exp(logRelativeWeight - relativeBoundKcal/RT) - 1 | logRelativeWeight > relativeBoundKcal/RT)',
        zero_tail_interpretation='undefined/no tail observed; not an estimated population zero',
        aggregation='per-state per-stage mean excess, then per-system maximum; distinct proposals not pooled',
        overlap_note='Discovery -initial checkpoints are subsets of their extended batch. All/prefinal counts sum recorded checkpoints and are not unique draw counts. Final batches are disjoint.',
        all=describe(artifacts), final=describe([r for r in artifacts if r['phase']=='final']),
        prefinal=describe([r for r in artifacts if r['phase']=='prefinal']),
        input_files=len(inventory), input_bytes=sum(r['bytes'] for r in inventory),
        output=str(out), errors=errors)
    table(out/'per_system.tsv', systems)
    table(out/'per_stage.tsv', artifacts)
    table(out/'exceeds_20.tsv', exceeds)
    table(out/'input_inventory.tsv', inventory)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    def fmt(v):
        return '无尾部观测' if v is None else f'{v:.6g}'
    lines = ['# Latest pair-only RB=1: observed conditional severity', '',
        'Source: '+str(ROOT), '',
        'Observed s0 = average(R/C - 1 among R > C). Per-system maxima over distinct state/stage summaries. No-tail batches are undefined, not evidence for zero population severity. These are descriptive empirical observations, not population upper bounds.', '',
        '| System | Final max | Prefinal max | Final batches with tail / all final | Final tail draws | All batches >20 |',
        '|---|---:|---:|---:|---:|---:|']
    for r in systems:
        lines.append(f"| {r['design']} | {fmt(r['final_max_s0'])} | {fmt(r['prefinal_max_s0'])} | {r['final_tail_batches']}/{r['final_batches']} | {r['final_tail_draws']} | {r['batches_over_20']} |")
    lines += ['', '## Every batch with observed s0 > 20', '',
              '| System | Stage | Observed s0 | Tail / N | State |', '|---|---|---:|---:|---|']
    for r in exceeds:
        lines.append(f"| {r['design']} | {r['stage']} | {r['observed_s0']:.9g} | {r['tail_count']}/{r['n']} | {r['state']} |")
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    print('SUMMARY='+json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
