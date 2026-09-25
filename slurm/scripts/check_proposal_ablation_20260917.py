"""End-to-end ablation gate; run only inside Slurm."""
import csv
import math
import os
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')

def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def keys(path):
    return {row['key']: row['value'] for row in table(path)}

root = Path(sys.argv[1]) / '1gwc' / 'seed_42'
baseline = Path('/usr/xtmp/lz280/packstar_forward_comparison_20260913/A12588114')
checks = []
for arm, old_arm in [('no-learning', None), ('pair-only', 'pair-only'),
                     ('decomposition-cost', 'current-selector'), ('budget-forward', 'budget-forward')]:
    run = root / arm
    result = table(run / '1gwc_packstar_pfunc.tsv')
    assert len(result) == 1, (arm, result)
    assert result[0]['status'] == 'Estimated', (arm, result)
    states = list((run / 'adaptive_frequency_severity').glob('state-*'))
    assert len(states) == 1, (arm, states)
    state = states[0]
    protocol = keys(state / 'frequency_severity_protocol.tsv')
    winner = keys(state / 'eta_selected.tsv')
    if arm == 'no-learning':
        assert protocol['proposalLearningEnabled'] == 'false'
        assert protocol['tripleEtaEnabled'] == 'false'
        assert winner['candidate'] == 'fixed-qm-no-learning'
        assert winner['proposalDpSweeps'] == '0'
        assert winner['tripleEtaSelectedPositionTriples'] == '0'
        log = (run / 'run.log').read_text()
        assert 'joint-pair-fit' not in log and 'triple-m2-fit' not in log
        sample_count = 0
        for path in state.glob('*_samples.tsv'):
            for sample in table(path):
                if 'etaKcal' not in sample:
                    continue
                assert float(sample['etaKcal']) == 0.0, (path, sample)
                assert float(sample['eMinKcal']) == float(sample['eProposalKcal']), (path, sample)
                sample_count += 1
        assert sample_count > 500, sample_count
        checks.append(f'no-learning: {sample_count} audited sample rows; eta=0, original proposal energy, zero corrected DP sweeps, no pair/triple fits')
    else:
        assert protocol['proposalLearningEnabled'] == 'true'
        old_run = baseline / old_arm / 'runs' / '1gwc_A12588114_T2'
        old_result = table(old_run / '1gwc_packstar_pfunc.tsv')[0]
        for key in ('epsilon', 'lower_log10', 'upper_log10'):
            assert math.isclose(float(result[0][key]), float(old_result[key]), rel_tol=0, abs_tol=1e-7), (arm, key, result, old_result)
        old_state = next((old_run / 'adaptive_frequency_severity').glob('state-*'))
        current_final = table(state / 'frequency_severity_final.tsv')[0]
        old_final = table(old_state / 'frequency_severity_final.tsv')[0]
        assert current_final['n'] == old_final['n'], (arm, current_final['n'], old_final['n'])
        assert winner['candidate'] == keys(old_state / 'eta_selected.tsv')['candidate']
        checks.append(f'{arm}: existing selected model/final N unchanged; interval agrees with frozen September 13 result to 1e-7')

(Path(sys.argv[1]) / 'gate_checks.txt').write_text('\n'.join(checks) + '\n')
print('\n'.join(checks))
