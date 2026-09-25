#!/bin/bash
#SBATCH --job-name=ps_ablation_gate
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:titan_v:4
#SBATCH --exclude=jerry3
#SBATCH --time=02:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${BUILD_ROOT:?}"
export CASE_OVERRIDE=1gwc SEED_OVERRIDE=42
export RESULT_ROOT_OVERRIDE=/usr/xtmp/lz280/packstar_proposal_ablation_20260917/gate_${SLURM_JOB_ID}
export PYTHONPYCACHEPREFIX="$RESULT_ROOT_OVERRIDE/pycache"
/bin/bash "$BUILD_ROOT/source/slurm/scripts/run_proposal_ablation_titan_20260917.slurm"
python3 "$BUILD_ROOT/source/slurm/scripts/check_proposal_ablation_20260917.py" "$RESULT_ROOT_OVERRIDE"
touch "$RESULT_ROOT_OVERRIDE/PASSED"
