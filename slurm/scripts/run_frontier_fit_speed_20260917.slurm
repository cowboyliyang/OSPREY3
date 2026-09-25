#!/bin/bash
#SBATCH --job-name=ps_frontier_fit_speed
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:a5000:8
#SBATCH --nodelist=fennario-02
#SBATCH --time=2-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${BUILD_ROOT:?}"
export SLURM_ARRAY_TASK_ID=0 SLURM_ARRAY_JOB_ID="$SLURM_JOB_ID"
RESULT=/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A${SLURM_JOB_ID}/runs/2xgy_A${SLURM_JOB_ID}_T0
mkdir -p "$RESULT"
# Freeze the small verification/entry scripts before executing. Project code
# and the complete 39-sequence runner already live in the validated build.
cp "${BASH_SOURCE[0]}" "$RESULT/speed_wrapper.slurm"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/check_frontier_fit_speed_20260917.py "$RESULT/check_speed.py"
/bin/bash "$BUILD_ROOT/source/slurm/scripts/run_frontier_fullnode_20260917.slurm"
export PYTHONPYCACHEPREFIX="$RESULT/pycache"
python3 "$RESULT/check_speed.py" \
    /usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A12623257/runs/2xgy_A12623257_T0 "$RESULT"
