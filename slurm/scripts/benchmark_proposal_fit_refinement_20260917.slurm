#!/bin/bash
#SBATCH --job-name=ps_fit_refine
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --nodelist=jerry6
#SBATCH --time=02:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${NEW_BUILD:?}"
BASELINE=/usr/xtmp/lz280/packstar_fit_speed_20260917/compare_12626605
REFINE_ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/refine_${SLURM_JOB_ID}
mkdir -p "$REFINE_ROOT"
export ARM_FILTER=budget-forward SEED_OVERRIDE=42 BUILD_ROOT="$NEW_BUILD"
export ALLOW_UNUSED_FOLD_STOPS=1 PYTHONPYCACHEPREFIX="$REFINE_ROOT/pycache"
for CASE_OVERRIDE in 1gwc 2xgy; do
    export CASE_OVERRIDE FIT_THREADS=16 RESULT_ROOT_OVERRIDE="$REFINE_ROOT/threads16"
    /bin/bash "$NEW_BUILD/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
    for BASELINE_MODE in legacy optimized; do
        python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_speed_20260917.py" \
            "$BASELINE/$BASELINE_MODE" "$RESULT_ROOT_OVERRIDE" "$CASE_OVERRIDE" \
            "$REFINE_ROOT/${CASE_OVERRIDE}_vs_${BASELINE_MODE}.json"
    done
done
export FIT_THREADS=32 RESULT_ROOT_OVERRIDE="$REFINE_ROOT/threads32"
/bin/bash "$NEW_BUILD/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_speed_20260917.py" \
    "$REFINE_ROOT/threads16" "$RESULT_ROOT_OVERRIDE" 2xgy "$REFINE_ROOT/2xgy_16_vs_32.json"
