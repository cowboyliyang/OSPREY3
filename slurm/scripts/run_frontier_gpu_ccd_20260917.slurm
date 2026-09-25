#!/bin/bash
#SBATCH --job-name=ps_frontier_gpu_ccd
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --array=0-5%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:a5000:8
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail

# Reuse the exact validated Java build and experiment protocol from the latest
# CPU CCD runs. Only the CCD backend changes; DP remains FP64 on eight GPUs.
export BUILD_ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626762
ROOT=/usr/xtmp/lz280/packstar_frontier_gpu_ccd_20260917/A${SLURM_ARRAY_JOB_ID:?}
TASK_ROOT="$ROOT/protocols/T${SLURM_ARRAY_TASK_ID:?}"
mkdir -p "$TASK_ROOT/tmp"
export TMPDIR="$TASK_ROOT/tmp"
export PYTHONPYCACHEPREFIX="$TASK_ROOT/pycache"
cp "${BASH_SOURCE[0]}" "$TASK_ROOT/submitted.slurm"

python3 - "$BUILD_ROOT" "$TASK_ROOT" <<'PY'
from pathlib import Path
import sys
build, task = map(Path, sys.argv[1:])
source = build/'source/slurm/scripts/run_frontier_fullnode_20260917.slurm'
text = source.read_text()
replacements = {
    '/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A':
        '/usr/xtmp/lz280/packstar_frontier_gpu_ccd_20260917/A',
    '-Dosprey.wmb.numGpus=0': '-Dosprey.wmb.numGpus=8',
    "printf 'ccdMode\\tcpu\\n'":
        "printf 'ccdMode\\tResidueCudaCCD\\n'\n"
        "    printf 'ccdPrecision\\tFP64-double-kernel\\n'\n"
        "    printf 'ccdGpuCount\\t8\\n'\n"
        "    printf 'ccdStreamsPerGpu\\t64\\n'",
}
for old, new in replacements.items():
    assert text.count(old) == 1, (old, text.count(old))
    text = text.replace(old, new)
(task/'runner.slurm').write_text(text)
(task/'protocol_changes.txt').write_text(
    'Baseline build: 12626762; CPU CCD baselines: 12626765 and 12627258.\n'
    'Changed osprey.wmb.numGpus from 0 to 8, selecting ResidueCudaCCD.\n'
    'Retained 64 streams/GPU, FP64 DP, all learning/sampling parameters and historical emat caches.\n'
    'Energy kernel and CCD coordinates/line search use double.\n')
PY

source /etc/profile.d/modules.sh 2>/dev/null || true
module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
KERNEL="$BUILD_ROOT/runtime/main-resources/gpuKernels/cuda/residueCcd.bin"
sha256sum "$KERNEL" \
    "$BUILD_ROOT/source/src/main/resources/gpuKernels/cuda/residueCcd.cu" \
    "$BUILD_ROOT/source/src/main/java/edu/duke/cs/osprey/gpu/cuda/kernels/ResidueCudaCCDMinimizer.java" \
    > "$TASK_ROOT/ccd_precision_inputs.sha256"
if command -v cuobjdump >/dev/null; then
    # This validated fatbin contains native sm_70/sm_86 code without PTX.
    cuobjdump --dump-sass "$KERNEL" > "$TASK_ROOT/residueCcd.sass.txt"
    if ! rg -m 8 '\bD(ADD|MUL|FMA)\b' "$TASK_ROOT/residueCcd.sass.txt" > "$TASK_ROOT/fp64_instructions.txt"; then
        echo 'No expected FP64 arithmetic instructions in CCD binary' >&2
        exit 6
    fi
fi
/bin/bash -n "$TASK_ROOT/runner.slurm"
/bin/bash "$TASK_ROOT/runner.slurm"
