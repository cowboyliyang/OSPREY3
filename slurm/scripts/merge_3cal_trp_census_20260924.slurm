#!/bin/bash
#SBATCH --job-name=ccd3cal_merge
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 - "${RESULT_ROOT:?}/shard_A${ARRAY_ID:?}" "${TARGET:?}" <<'PY'
import sys
from pathlib import Path
root=Path(sys.argv[1]);target=sys.argv[2]
for task in range(931):
    assert (root/f'T{task}'/'target.txt').read_text().strip()==target, f'target mismatch {task}'
PY
python3 "${BUILD_ROOT:?}/source/merge_million_census_20260918.py" "${RESULT_ROOT:?}" "${ARRAY_ID:?}" "$BUILD_ROOT/config.tsv" "${BLOCK:-32768}"
python3 - "$RESULT_ROOT/merged_A$ARRAY_ID" "$TARGET" <<'PY'
import json,math,sys
from pathlib import Path
out=Path(sys.argv[1]); report=json.loads((out/'reference.json').read_text())
z=report['logZ']/math.log(10)
target=sys.argv[2]
intervals={'PACK_rb1': [70.950631,71.311520]} if target=='sparse' else {'historical_MARK': [71.435162,71.925444]}
comparison={'log10_Z':z,'sequence':'TRP','state':'protein','target':target,'intervals':intervals,
            'inside':{k:lo<=z<=hi for k,(lo,hi) in intervals.items()}}
(out/'comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
with (out/'reference.md').open('a') as f:
    f.write(f'\n3CAL TRP / protein, target={target}: log10(Z) = {z:.14g}.\n\n')
    for k,(lo,hi) in intervals.items():
        f.write(f'{k}: [{lo}, {hi}], contains exhaustive CCD value: {lo<=z<=hi}.\n\n')
print(json.dumps(comparison,indent=2))
PY
