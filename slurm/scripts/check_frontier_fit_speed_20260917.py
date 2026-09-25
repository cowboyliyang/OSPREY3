"""Exact full-workload check; only timings and unused fold-stop metadata vary."""
import csv
import json
import os
import re
import sys
from pathlib import Path

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run through Slurm.')

def rows(path, delimiter='\t'):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))

old, new = map(Path, sys.argv[1:3])
before, after = [rows(p / '2xgy_packstar.csv', ',') for p in (old, new)]
assert len(before) == len(after) == 39
for left, right in zip(before, after):
    assert left.keys() == right.keys()
    for key in left:
        if key != 'total_time_s':
            assert left[key] == right[key], (left['sequence'], key, left[key], right[key])

roots = [p / 'adaptive_frequency_severity' for p in (old, new)]
files = [{p.relative_to(root) for p in root.glob('state-*/*.tsv')} for root in roots]
assert files[0] == files[1], 'state/artifact sets changed'
checked = 0
diagnostic_changes = []
for name in sorted(files[0]):
    if name.name == 'frequency_severity_protocol.tsv':
        continue
    paths = [root / name for root in roots]
    if paths[0].read_bytes() != paths[1].read_bytes():
        assert name.name.endswith('eta_candidate_model_path.tsv'), str(name)
        tables = [rows(p) for p in paths]
        assert len(tables[0]) == len(tables[1]), str(name)
        for left, right in zip(*tables):
            assert left.keys() == right.keys(), str(name)
            for key in left:
                if left[key] == right[key]:
                    continue
                assert key in ('fitFold0Stop', 'fitFold1Stop') and right[key] == 'required-prefix-limit', (str(name), key)
                diagnostic_changes.append(dict(artifact=str(name), candidate=left['candidate'], column=key))
    checked += 1
walls = [float(re.search(r'elapsed=([\d.]+)', (p / 'wall.time').read_text())[1]) for p in (old, new)]
manifests = [dict(line.split('\t', 1) for line in (p / 'run_manifest.tsv').read_text().splitlines()) for p in (old, new)]
for key in ('designId', 'node', 'precision', 'cpuThreads', 'gpuName', 'branchDpGpuCount', 'scope', 'targetEpsilon'):
    assert manifests[0][key] == manifests[1][key], key
report = dict(case='2xgy', scope='39-sequences-P-L-PL', same_node=True,
              exact_scored_models_samples_intervals=True, checked_artifacts=checked,
              changed_unused_fold_diagnostics=diagnostic_changes,
              old_wall_s=walls[0], new_wall_s=walls[1], speedup=walls[0]/walls[1],
              runs=[str(old), str(new)])
(new / 'full_workload_speed_check.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'changed_unused_fold_diagnostics'}, indent=2))
