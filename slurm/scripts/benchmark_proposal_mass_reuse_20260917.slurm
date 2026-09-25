#!/bin/bash
#SBATCH --job-name=ps_mass_reuse_gate
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --exclude=jerry3
#SBATCH --time=03:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${NEW_BUILD:?}"
BASE_BUILD=/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626676
ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/mass_reuse_${SLURM_JOB_ID}
mkdir -p "$ROOT"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
export ARM_FILTER=budget-forward SEED_OVERRIDE=42 FIT_THREADS=16
export PYTHONPYCACHEPREFIX="$ROOT/pycache"
RUNNER="$NEW_BUILD/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
# Matched hardware, frozen baseline, identical seed/protocol. Alternate order.
index=0
for CASE_OVERRIDE in 2xgy 1gwc 4wyu 2xgy; do
    export CASE_OVERRIDE SEQ_OVERRIDE=0
    if [[ "$index" == 3 ]]; then export SEQ_OVERRIDE=36; fi
    pair="$ROOT/${CASE_OVERRIDE}_seq${SEQ_OVERRIDE}"
    modes=(base reused)
    if (( index % 2 )); then modes=(reused base); fi
    for mode in "${modes[@]}"; do
        export RESULT_ROOT_OVERRIDE="$pair/$mode"
        if [[ "$mode" == base ]]; then export BUILD_ROOT="$BASE_BUILD"; else export BUILD_ROOT="$NEW_BUILD"; fi
        /bin/bash "$RUNNER"
    done
    python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_numerics_20260917.py" \
        "$pair/base" "$pair/reused" "$CASE_OVERRIDE" "$pair/check.json"
    index=$((index + 1))
done
touch "$ROOT/PASSED"
