#!/bin/bash
#SBATCH --job-name=ps_triple_profile
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --exclude=jerry3
#SBATCH --time=01:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export BUILD_ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626676
ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/profile_${SLURM_JOB_ID}
mkdir -p "$ROOT/tmp"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
export ARM_FILTER=budget-forward SEED_OVERRIDE=42 FIT_THREADS=16 CASE_OVERRIDE=2xgy SEQ_OVERRIDE=36
export RESULT_ROOT_OVERRIDE="$ROOT/run"
export JAVA_TOOL_OPTIONS="-XX:StartFlightRecording=filename=$ROOT/cpu.jfr,settings=profile,dumponexit=true"
/bin/bash "$BUILD_ROOT/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
unset JAVA_TOOL_OPTIONS
/home/users/lz280/java/jdk-17.0.2+8/bin/jfr print --json --events jdk.ExecutionSample "$ROOT/cpu.jfr" > "$ROOT/samples.json"
python3 - "$ROOT" <<'PY'
import collections, json, sys
from pathlib import Path
root = Path(sys.argv[1])
events = json.loads((root/'samples.json').read_text())['recording']['events']
top, inclusive = collections.Counter(), collections.Counter()
total = 0
for event in events:
    trace = event['values'].get('stackTrace')
    if not trace: continue
    frames = trace['frames']
    names = [f['method']['type']['name'].replace('/', '.') + '.' + f['method']['name'] for f in frames]
    if not names: continue
    total += 1
    top[names[0]] += 1
    inclusive.update(set(names))
lines = [f'Execution samples: {total}', 'Self samples:']
lines += [f'{n:6d} {name}' for name, n in top.most_common(35)]
lines += ['Inclusive samples (overlap; not wall-time fractions):']
lines += [f'{n:6d} {name}' for name, n in inclusive.most_common(45)]
(root/'profile.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
PY
