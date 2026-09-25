"""Compare deterministic scientific outputs and report measured process time."""
import csv
import json
import os
import re
import sys
from pathlib import Path
if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')

def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def keys(path):
    return {r['key']: r['value'] for r in table(path)}

old_root, new_root, case, output = sys.argv[1:]
runs = [Path(p) / case / 'seed_42' / 'budget-forward' for p in (old_root, new_root)]
results = [table(p / f'{case}_packstar_pfunc.tsv')[0] for p in runs]
for key in ('status', 'epsilon', 'lower_log10', 'upper_log10', 'sequence', 'num_confs'):
    assert results[0][key] == results[1][key], (key, results)
states = [next((p / 'adaptive_frequency_severity').glob('state-*')) for p in runs]
assert {p.name for p in states[0].glob('*.tsv')} == {p.name for p in states[1].glob('*.tsv')}
# All statistical/model/sample artifacts must match byte for byte. Protocol
# adds the performance-thread setting in the optimized version only.
checked = []
changed_diagnostics = []
for old in sorted(states[0].glob('*.tsv')):
    if old.name == 'frequency_severity_protocol.tsv':
        continue
    new = states[1] / old.name
    if old.read_bytes() != new.read_bytes():
        assert os.environ.get('ALLOW_UNUSED_FOLD_STOPS') == '1' and old.name.endswith('eta_candidate_model_path.tsv'), f'changed scientific artifact: {old.name}'
        before, after = table(old), table(new)
        assert len(before) == len(after), old.name
        for left, right in zip(before, after):
            assert left.keys() == right.keys(), old.name
            for key in left:
                if left[key] == right[key]:
                    continue
                assert key in ('fitFold0Stop', 'fitFold1Stop') and right[key] == 'required-prefix-limit', (old.name, key, left[key], right[key])
                changed_diagnostics.append(dict(artifact=old.name, candidate=left['candidate'], column=key, before=left[key], after=right[key]))
    checked.append(old.name)
walls = []
nodes = []
for run in runs:
    walls.append(float(re.search(r'elapsed=([\d.]+)', (run / 'wall.time').read_text())[1]))
    manifest = dict(line.split('\t', 1) for line in (run / 'run_manifest.tsv').read_text().splitlines() if '\t' in line)
    nodes.append(manifest['node'])
path_ms = sum(float(x) for x in re.findall(r'pathMs=([\d.Ee+-]+)', (runs[1] / 'run.log').read_text()))
report = dict(case=case, status=results[1]['status'], checked_artifacts=len(checked),
              exact_scientific_artifact_agreement=not changed_diagnostics,
              exact_scored_models_samples_intervals=True, changed_unused_fold_diagnostics=changed_diagnostics,
              old_wall_s=walls[0], new_wall_s=walls[1],
              observed_wall_ratio=walls[0]/walls[1], old_node=nodes[0], new_node=nodes[1],
              same_node=nodes[0] == nodes[1], optimized_path_s=path_ms/1000,
              runs=[str(p) for p in runs])
Path(output).write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
