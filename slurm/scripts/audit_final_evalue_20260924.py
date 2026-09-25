"""Replay frozen rb=1,2,4 severity diagnostics; run only through Slurm."""
import csv
import json
import math
import os
from pathlib import Path

RB1 = Path('/usr/xtmp/lz280/packstar_rb1_s0_audit_20260923_12690121/per_stage.tsv')
RB24 = Path('/usr/xtmp/lz280/packstar_rb24_final_s0_audit_20260924_12690389/all_final.tsv')
RT = 0.0019891 * 298.15
LAMBDAS = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, .1, .3, .5, .8)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def keys(path):
    return {r['key']: r['value'] for r in rows(path)}


def logsum(values):
    high = max(values)
    return high + math.log(math.fsum(math.exp(v - high) for v in values))


def log_evalue(log_excesses, cap):
    if len(log_excesses) <= 1:
        return math.nan
    ratios = [v - math.log(cap) for v in log_excesses]
    components = [logsum(ratios) - math.log(len(ratios))]
    for lam in LAMBDAS:
        components.append(math.fsum(logsum((math.log1p(-lam), math.log(lam) + r)) for r in ratios))
    return logsum(components) - math.log(len(components))


def write_table(path, records):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def main():
    job = os.environ['SLURM_JOB_ID']
    out = Path('/usr/xtmp/lz280') / ('packstar_final_evalue_audit_20260924_' + job)
    out.mkdir(exist_ok=False)
    records, replay = [], []
    inputs = [dict(r, rb=1) for r in rows(RB1) if r['stage'] == 'final'] + rows(RB24)
    for source in inputs:
        path = Path(source['path'])
        stage = rows(path)[0]
        protocol = keys(path.with_name('frequency_severity_protocol.tsv'))
        monitor = rows(path.with_name('frequency_severity_monitor.tsv'))[0]
        cap = float(protocol['severityCapS0'])
        log_clip = float(protocol['relativeBoundKcal']) / RT
        alpha = float(protocol['severityTestAlpha'])
        assert cap == 20 and alpha == .05
        tail_count = int(stage['tailCount'])
        log_excesses, tail_weights, tail_confs = [], [], set()
        if tail_count:
            samples = rows(path.with_name('frequency_severity_final_samples.tsv'))
            assert len(samples) == int(stage['n'])
            for sample in samples:
                log_weight = float(sample['logRelativeWeight'])
                relative = log_weight - log_clip
                assert (relative > 0) == (sample['tailEvent'] == 'true')
                if relative > 0:
                    log_excesses.append(relative + math.log1p(-math.exp(-relative))
                                        if relative > 40 else math.log(math.expm1(relative)))
                    tail_weights.append(log_weight)
                    tail_confs.add(sample['conf'])
        assert len(log_excesses) == tail_count
        mean = float(stage['empiricalConditionalSeverity'])
        if tail_count:
            assert math.isclose(logsum(log_excesses) - math.log(tail_count), math.log(mean), abs_tol=1e-9)
        log_e = log_evalue(log_excesses, cap)
        monitor_log_e = float(monitor['severityTestLogE'])
        final_reject = log_e >= -math.log(alpha)
        rb = int(source['rb'])
        case_id = f"rb{rb}:{source['design']}:{path.parent.name}"
        record = dict(rb=rb, design=source['design'], case_id=case_id,
                      sequence=source.get('sequence', ''), role=source.get('role', ''),
                      result_matches=source.get('final_result_matches', ''),
                      result_status=source.get('result_status', ''),
                      n=int(stage['n']), tail_count=tail_count,
                      tail_unique_confs=len(tail_confs), empirical_severity=mean,
                      individual_severity_over_cap=sum(v > math.log(cap) for v in log_excesses),
                      monitor_tail=int(monitor['tailCount']), monitor_log_e=monitor_log_e,
                      final_log_e=log_e,
                      final_evalue=math.nan if math.isnan(log_e) else math.exp(log_e) if log_e < 709 else math.inf,
                      test_alpha=alpha,
                      old_monitor_rejected=monitor['severityTestRejected'],
                      new_final_rejected=final_reject, path=str(path))
        records.append(record)
        # Bulk values do not enter this e-value. Pad to satisfy the Java API's
        # minimum batch length while preserving all tail multiplicities.
        compact_weights = tail_weights + [log_clip - 1] * max(0, 2 - tail_count)
        replay.append(dict(case_id=case_id, log_clip=log_clip, cap=cap, alpha=alpha,
                           tail_count=tail_count, expected_log_e=log_e,
                           expected_rejected=str(final_reject).lower(),
                           log_weights=','.join(map(repr, compact_weights))))
    write_table(out / 'all_final.tsv', records)
    write_table(out / 'java_replay.tsv', replay)
    notable = sorted((r for r in records if r['empirical_severity'] > 20 or r['new_final_rejected']),
                     key=lambda r: (r['rb'], -r['empirical_severity']))
    write_table(out / 'notable.tsv', notable)
    summary = {str(rb): dict(final_batches=sum(r['rb'] == rb for r in records),
                            final_with_tail=sum(r['rb'] == rb and r['tail_count'] > 0 for r in records),
                            final_rejected=sum(r['rb'] == rb and r['new_final_rejected'] for r in records))
               for rb in (1, 2, 4)}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    lines = ['# Replay of a single fixed final e-value test', '',
             'The existing 11-component formula is unchanged. One final test at alpha=0.05 rejects at e>=20; no separate monitor. NaN means fewer than two tail samples (inconclusive). No new CCD evaluations.', '',
             '| rb | System | Sequence/state match | Severity | Tail/N | Tail values >20 | Monitor tail | Final e | Final reject |',
             '|---|---|---|---:|---:|---:|---:|---:|---|']
    for r in notable:
        label = r['result_matches'] or (r['sequence'] + ' / ' + r['role'])
        lines.append(f"| {r['rb']} | {r['design']} | {label} | {r['empirical_severity']:.6g} | {r['tail_count']}/{r['n']} | {r['individual_severity_over_cap']} | {r['monitor_tail']} | {r['final_evalue']:.6g} | {r['new_final_rejected']} |")
    (out / 'report.md').write_text('\n'.join(lines) + '\n')
    print((out / 'report.md').read_text(), flush=True)
    print(json.dumps(summary), flush=True)
    print('OUTPUT=' + str(out), flush=True)


if __name__ == '__main__':
    main()
