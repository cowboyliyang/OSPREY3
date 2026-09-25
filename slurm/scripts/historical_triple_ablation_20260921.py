"""Current-code ablation on the audited August triple-development hard cases.

Reuse the matched experiment runner; replace only historical selection and
input preparation. Every executable entry point remains Slurm-only.
"""
import csv
import json
from pathlib import Path
import shutil

import hardcase_ablation_20260921 as experiment


SPEC = Path(__file__).with_name('historical_triple_cases_20260921.json')


def metadata(path):
    return {row.get('key', row.get('field')): row['value']
            for row in experiment.table(path)}


def historical_evidence(case):
    pair_run = Path(case['historical_pair_run'])
    triple_run = Path(case['historical_triple_run'])
    pair_state = pair_run / case['historical_pair_state_dir']
    triple_state = triple_run / case['historical_triple_state_dir']
    assert pair_state.name == triple_state.name
    result_name = f'{case["design"]}_{case["state"]}_seq{case["seq_index"]}_packstar_pfunc.tsv'
    paths = {
        'pair_result': pair_run / result_name,
        'triple_result': triple_run / result_name,
        'pair_manifest': pair_run / 'run_manifest.tsv',
        'triple_manifest': triple_run / 'run_manifest.tsv',
        'pair_protocol': pair_state / 'frequency_severity_protocol.tsv',
        'triple_protocol': triple_state / 'frequency_severity_protocol.tsv',
        'triple_final': triple_state / 'frequency_severity_final_interval.tsv',
    }
    evidence = {'case': case['id'], 'files': {
        label: {'path': str(path), 'sha256': experiment.sha(path)}
        for label, path in paths.items()}}
    for label in ('pair_result', 'triple_result'):
        rows = experiment.table(paths[label])
        assert len(rows) == 1, (case, label, rows)
        row = rows[0]
        assert row['design_id'] == case['design']
        assert row['state'].lower() == case['state']
        assert int(row['seq_index']) == case['seq_index']
        assert row['sequence'].strip() == case['sequence']
        assert row['status'] == ('Aborted' if label == 'pair_result' else 'Estimated')
        if label == 'triple_result':
            assert float(row['epsilon']) <= 0.683
        evidence[label] = row
    for label in ('pair_manifest', 'triple_manifest'):
        values = metadata(paths[label])
        assert values['designId'] == case['design']
        assert values['pfuncState'] == case['state']
        assert int(values['pfuncSequenceIndex']) == case['seq_index']
        for key, expected in [('branchResidualBudget', 1), ('randomSeed', 42),
                              ('trainSamples', 500), ('finalMaxSamples', 4000),
                              ('targetEpsilon', 0.683), ('confidenceDelta', 0.05)]:
            assert float(values[key]) == expected, (case, label, key, values[key])
        evidence[label] = values
    pair = metadata(paths['pair_protocol'])
    triple = metadata(paths['triple_protocol'])
    # The legacy schema predates the triple switch; eta-only is its explicit
    # algorithm declaration, not an inference from an Aborted result.
    assert 'V2.1-eta-only' in pair['algorithm'], pair
    assert pair.get('tripleEtaEnabled', 'false') == 'false'
    assert triple['tripleEtaEnabled'] == 'true'
    assert pair['randomStreamIdentityHash'] == triple['randomStreamIdentityHash']
    assert pair_state.name == 'state-' + pair['randomStreamIdentityHash']
    for protocol in (pair, triple):
        for key, expected in [('trainSamples', 500), ('maxFinalSamples', 4000),
                              ('targetEpsilon', 0.683), ('confidenceDelta', 0.05),
                              ('relativeBoundKcal', 1), ('severityCapS0', 20)]:
            assert float(protocol[key]) == expected, (case, key, protocol[key])
    final = metadata(paths['triple_final'])
    for key in ('selectedTripleEtaActive', 'targetReached', 'certificateValid'):
        assert final[key] == 'true', (case, key, final[key])
    assert int(final['selectedTripleEtaPositionTriples']) == 3
    assert float(final['epsilon']) <= 0.683
    evidence.update(pair_protocol=pair, triple_protocol=triple, triple_final=final)
    failure_path = pair_state / 'frequency_severity_failure.tsv'
    if failure_path.exists():
        evidence['files']['pair_failure'] = {'path': str(failure_path),
                                             'sha256': experiment.sha(failure_path)}
        evidence['pair_failure'] = metadata(failure_path)
    return evidence


def prepare(build):
    spec = json.loads(SPEC.read_text())
    assert spec['seeds'] == [42, 43, 44]
    assert spec['arms'] == ['no-learning', 'pair-only', 'budget-forward']
    assert all(case['residual_budget'] == '1' for case in spec['cases'])
    historical = [historical_evidence(case) for case in spec['cases']]
    inputs = build / 'inputs'
    inputs.mkdir()
    csv_path = Path('/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv')
    designs = {case['design'] for case in spec['cases']}
    rows = {row[0]: row for row in csv.reader(csv_path.read_text().splitlines())
            if row and row[0] in designs}
    assert set(rows) == designs
    input_hashes = {}
    sequence_maps = {}
    for design, row in rows.items():
        dest = inputs / design
        dest.mkdir()
        pdb = Path('/usr/xtmp/lz280/dance_bench/pdbs_prepped') / row[1] / (row[1] + '.min.reduce.renum.pdb')
        sources = {dest / 'input.pdb': pdb}
        for state in ('complex', 'protein', 'ligand'):
            for old, new in (('rigid', 'rigid'), ('min', 'minimizing')):
                source = Path('/usr/xtmp/lz280/bench_comparison/results/emat_cache') / design / f'export.{state}.{old}.dat'
                sources[dest / f'packstar.{state}.{new}.dat'] = source
        for target, source in sources.items():
            shutil.copy2(source, target)
            digest = experiment.sha(source)
            assert experiment.sha(target) == digest
            input_hashes[str(target.relative_to(inputs))] = {'source': str(source), 'sha256': digest}
        (dest / 'design_row.csv').write_text(','.join(row) + '\n')
        props = {'osprey.bench.method': 'sequence_dump', 'osprey.bench.outputDir': dest,
                 'osprey.bench.designId': design, 'osprey.bench.pdbPath': dest / 'input.pdb',
                 'osprey.bench.mutable': row[5], 'osprey.bench.flexible': row[6],
                 'osprey.bench.numCPUs': 1, 'osprey.sequenceDump.maxMut': 1,
                 'osprey.sequenceDump.output': dest / 'sequences.tsv',
                 'java.io.tmpdir': build / 'tmp'}
        experiment.checked_run(experiment.java_command(build, props, '8g'),
                               dest / 'sequence_dump.log', cwd=build / 'source')
        sequence_maps[design] = {int(row['seq_index']): row['sequence'].strip()
                                 for row in experiment.table(dest / 'sequences.tsv')}
    for case in spec['cases']:
        assert sequence_maps[case['design']][case['seq_index']] == case['sequence']
        case['mutable'], case['flexible'] = rows[case['design']][5:7]
    experiment.dump(build / 'cases.json', spec)
    experiment.dump(build / 'input_manifest.json', input_hashes)
    experiment.dump(build / 'historical_evidence.json', {
        'selection': spec['selection'],
        'cases': historical})
    print(f'PREPARED: {len(historical)} explicit eta-only failures and FINAL K=3 successes; current sequence identities verified.', flush=True)


original_summarize = experiment.summarize


def summarize(build, root):
    original_summarize(build, root)
    report = root / 'comparison.md'
    text = report.read_text().replace(
        'Historical failed-then-successful-with-final-triples observations selected the cases; they are not matched controls.',
        'The three August development cases have explicit eta-only failure protocols and later FINAL K=3 successes. Historical versions also differ in adaptation and overlap gates; only the current frozen-code arms are matched ablation controls.')
    report.write_text(text)


if __name__ == '__main__':
    experiment.prepare = prepare
    experiment.summarize = summarize
    experiment.main()
