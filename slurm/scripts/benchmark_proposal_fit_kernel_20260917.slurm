#!/bin/bash
#SBATCH --job-name=ps_fit_kernel
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --nodelist=jerry6
#SBATCH --time=03:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${NEW_BUILD:?}"
BASE_BUILD=/usr/xtmp/lz280/packstar_fit_speed_20260917/build_12626642
BASE_RESULTS=/usr/xtmp/lz280/packstar_fit_speed_20260917/refine_12626643/threads16
KERNEL_ROOT=/usr/xtmp/lz280/packstar_fit_speed_20260917/kernel_${SLURM_JOB_ID}
mkdir -p "$KERNEL_ROOT"
export ARM_FILTER=budget-forward SEED_OVERRIDE=42 FIT_THREADS=16
export PYTHONPYCACHEPREFIX="$KERNEL_ROOT/pycache"
RUNNER="$NEW_BUILD/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
for CASE_OVERRIDE in 1gwc 2xgy 4wyu; do
    export CASE_OVERRIDE SEQ_OVERRIDE=0
    if [[ "$CASE_OVERRIDE" == 4wyu ]]; then
        export BUILD_ROOT="$BASE_BUILD" RESULT_ROOT_OVERRIDE="$KERNEL_ROOT/base"
        /bin/bash "$RUNNER"
        OLD_RESULTS="$RESULT_ROOT_OVERRIDE"
    else
        OLD_RESULTS="$BASE_RESULTS"
    fi
    export BUILD_ROOT="$NEW_BUILD" RESULT_ROOT_OVERRIDE="$KERNEL_ROOT/shifted"
    /bin/bash "$RUNNER"
    python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_numerics_20260917.py" \
        "$OLD_RESULTS" "$RESULT_ROOT_OVERRIDE" "$CASE_OVERRIDE" "$KERNEL_ROOT/${CASE_OVERRIDE}_check.json"
done
export CASE_OVERRIDE=2xgy SEQ_OVERRIDE=36
for MODE in base shifted; do
    export RESULT_ROOT_OVERRIDE="$KERNEL_ROOT/mutant36_$MODE"
    if [[ "$MODE" == base ]]; then export BUILD_ROOT="$BASE_BUILD"; else export BUILD_ROOT="$NEW_BUILD"; fi
    /bin/bash "$RUNNER"
done
python3 "$NEW_BUILD/source/slurm/scripts/check_proposal_fit_numerics_20260917.py" \
    "$KERNEL_ROOT/mutant36_base" "$KERNEL_ROOT/mutant36_shifted" 2xgy "$KERNEL_ROOT/2xgy_mutant36_check.json"
