#!/bin/bash
#SBATCH --job-name=ps_old_triple_fullnode
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:a5000:8
#SBATCH --nodelist=fennario-01
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export BUILD_ROOT=/usr/xtmp/lz280/packstar_rescue_20260905/current_build
export SLURM_ARRAY_TASK_ID=0 SLURM_ARRAY_JOB_ID="$SLURM_JOB_ID"
ROOT=/usr/xtmp/lz280/packstar_historical_triple_20260917/A${SLURM_JOB_ID}
mkdir -p "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
# Generate a new runner; never modify either frozen build or its scripts.
# Use the exact September 5 runtime, with today's full-node resource controls.
python3 - "$ROOT" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
source = Path('/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626676/source/slurm/scripts/run_frontier_fullnode_20260917.slurm')
text = source.read_text()
replacements = {
    'REPO="$BUILD_ROOT/source"': 'REPO="/home/users/lz280/OSPREY3-fresh-packstar"',
    'test -f "$BUILD_ROOT/READY"': 'test -s "$BUILD_ROOT/test_classpath.txt"',
    '/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A${ARRAY_JOB_ID}': '/usr/xtmp/lz280/packstar_historical_triple_20260917/A${ARRAY_JOB_ID}',
    'GIT_HEAD="$(cat "$BUILD_ROOT/git_head.txt")"': 'GIT_HEAD="1456ab8e5b74ddd8e74e757d073afbb20b07eb10"',
    'budget-forward': 'historical-sep05-conditional-moment',
}
for old, new in replacements.items():
    assert old in text, old
    text = text.replace(old, new)
# These controls were introduced after September 5. Use the old build's
# original selection logic, not a new build with an approximate legacy flag.
text = '\n'.join(line for line in text.splitlines()
                 if not any(key in line for key in (
                     '-Dpackstar.pac.frequencySeverity.tripleEtaSelectionStrategy=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxWork=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxMStates=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxTableBytes=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxTotalTableBytes=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxHostBytes=',
                     '-Dpackstar.pac.frequencySeverity.tripleEtaMaxFileBytes='))) + '\n'
(root/'runner.slurm').write_text(text)
PY
# Audit the actual frozen bytecode/resources and input classpath on the node.
rg --files "$BUILD_ROOT/gradle_build/osprey-root/classes" "$BUILD_ROOT/gradle_build/osprey-root/resources" -0 | sort -z | xargs -0 sha256sum > "$ROOT/runtime.sha256"
cp "$BUILD_ROOT/test_classpath.txt" "$ROOT/historical_classpath.txt"
/bin/bash -n "$ROOT/runner.slurm"
/bin/bash "$ROOT/runner.slurm"
