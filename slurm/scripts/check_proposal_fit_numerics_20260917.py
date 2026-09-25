"""FP64-reduction gate: exact choices/counts/draws, tightly bounded floats."""
import csv
import json
import math
import os
import re
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')

def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

old_root, new_root, case, output = sys.argv[1:]
runs = [Path(p) / case / 'seed_42' / 'budget-forward' for p in (old_root, new_root)]
states = [next((p / 'adaptive_frequency_severity').glob('state-*')) for p in runs]
assert states[0].name == states[1].name
assert {p.name for p in states[0].glob('*.tsv')} == {p.name for p in states[1].glob('*.tsv')}
numeric_changes = 0
hash_changes = 0
max_normalized_error = 0.0
checked = []
accepted_candidate_differences = []
numeric_failures = []

def compare(left, right, label):
    global numeric_changes, hash_changes, max_normalized_error
    assert len(left) == len(right), (label, 'row count')
    for a, b in zip(left, right):
        assert a.keys() == b.keys(), label
        for key in a:
            if a[key] == b[key] or key == 'elapsed_s':
                continue
            logical_key = a.get('key', key) if key == 'value' else key
            if logical_key in ('signature', 'signatureSha256', 'tripleEtaSignatureSha256'):
                assert re.fullmatch('[0-9a-f]{64}', a[key]) and re.fullmatch('[0-9a-f]{64}', b[key]), (label, logical_key)
                hash_changes += 1
                continue
            # Integer counts, K, sample/sequence indices and booleans remain
            # exact. The user accepted small unselected-candidate diagnostic
            # drift on September 17 after reviewing the 6.85e-8 discrepancy.
            assert not re.fullmatch(r'[+-]?\d+', a[key]), (label, logical_key, a[key], b[key])
            try:
                x, y = float(a[key]), float(b[key])
            except ValueError:
                raise AssertionError((label, logical_key, a[key], b[key]))
            assert math.isfinite(x) and math.isfinite(y), (label, logical_key)
            originally_close = math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-8)
            unselected_candidate = 'eta_candidate_' in label and a.get('selected') == 'false'
            accepted_diagnostic = unselected_candidate and math.isclose(x, y, rel_tol=1e-7, abs_tol=1e-7)
            if not originally_close:
                detail = dict(artifact=label, candidate=a.get('candidate'), column=logical_key,
                              before=x, after=y, absolute_difference=abs(x-y))
                if accepted_diagnostic:
                    accepted_candidate_differences.append(detail)
                else:
                    numeric_failures.append(detail)
            numeric_changes += 1
            max_normalized_error = max(max_normalized_error, abs(x-y) / max(1, abs(x), abs(y)))

compare(*[rows(p / f'{case}_packstar_pfunc.tsv') for p in runs], 'pfunc')
for old in sorted(states[0].glob('*.tsv')):
    if old.name == 'frequency_severity_protocol.tsv':
        continue
    new = states[1] / old.name
    if old.read_bytes() != new.read_bytes():
        compare(rows(old), rows(new), old.name)
    checked.append(old.name)
walls = [float(re.search(r'elapsed=([\d.]+)', (p / 'wall.time').read_text())[1]) for p in runs]
paths = [sum(float(x) for x in re.findall(r'pathMs=([\d.Ee+-]+)', (p / 'run.log').read_text()))/1000 for p in runs]
nodes = [dict(line.split('\t', 1) for line in (p/'run_manifest.tsv').read_text().splitlines())['node'] for p in runs]
report = dict(case=case, exact_choices_counts_draws=True, checked_artifacts=len(checked),
              absolute_tolerance=1e-8, relative_tolerance=1e-10, numeric_changes=numeric_changes,
              acceptance='user-approved-small-unselected-candidate-diagnostic-drift-2026-09-17',
              unselected_candidate_absolute_tolerance=1e-7, unselected_candidate_relative_tolerance=1e-7,
              accepted_candidate_differences=accepted_candidate_differences, numeric_failures=numeric_failures,
              passed=not numeric_failures,
              changed_float_hashes=hash_changes, max_normalized_error=max_normalized_error,
              old_wall_s=walls[0], new_wall_s=walls[1], wall_speedup=walls[0]/walls[1],
              old_path_s=paths[0], new_path_s=paths[1], path_speedup=paths[0]/paths[1],
              same_node=nodes[0] == nodes[1], nodes=nodes, runs=[str(p) for p in runs])
Path(output).write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
if numeric_failures:
    raise SystemExit('Unaccepted numeric differences remain; see report.')
