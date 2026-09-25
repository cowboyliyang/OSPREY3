#!/bin/bash
#SBATCH --job-name=ps_fit_speed
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --exclude=jerry3
#SBATCH --time=12:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${NEW_BUILD:?}"
BASE_BUILD=/usr/xtmp/lz280/packstar_proposal_ablation_20260917/build_12624717
SPEED_ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/compare_${SLURM_JOB_ID}
mkdir -p "$SPEED_ROOT"
export ARM_FILTER=budget-forward SEED_OVERRIDE=42
export PYTHONPYCACHEPREFIX="$SPEED_ROOT/pycache"
RUNNER="$NEW_BUILD/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
# A first new-build run compares against the just-completed frozen real-system
# gate. Then a second case runs old/new sequentially on this same node.
export CASE_OVERRIDE=1gwc BUILD_ROOT="$NEW_BUILD" FIT_THREADS=16
export RESULT_ROOT_OVERRIDE="$SPEED_ROOT/optimized"
/bin/bash "$RUNNER"
python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_speed_20260917.py" \
    /usr/xtmp/lz280/packstar_proposal_ablation_20260917/gate_12624718 \
    "$SPEED_ROOT/optimized" 1gwc "$SPEED_ROOT/first_check.json"
export BUILD_ROOT="$BASE_BUILD" RESULT_ROOT_OVERRIDE="$SPEED_ROOT/legacy"
/bin/bash "$RUNNER"
python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_speed_20260917.py" \
    "$SPEED_ROOT/legacy" "$SPEED_ROOT/optimized" 1gwc "$SPEED_ROOT/1gwc_check.json"
for CASE_OVERRIDE in 2xgy; do
    export CASE_OVERRIDE
    for MODE in legacy optimized; do
        export RESULT_ROOT_OVERRIDE="$SPEED_ROOT/$MODE"
        if [[ "$MODE" == legacy ]]; then export BUILD_ROOT="$BASE_BUILD"; else export BUILD_ROOT="$NEW_BUILD"; fi
        /bin/bash "$RUNNER"
    done
    python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_speed_20260917.py" \
        "$SPEED_ROOT/legacy" "$SPEED_ROOT/optimized" "$CASE_OVERRIDE" "$SPEED_ROOT/${CASE_OVERRIDE}_check.json"
done
