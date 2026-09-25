#!/bin/bash
#SBATCH --job-name=pc-experimental-prep
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
RUN_ROOT="/usr/xtmp/lz280/packstar_functional_validation/experimental_${SLURM_JOB_ID}"
mkdir -p "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export PYTHONPYCACHEPREFIX="$RUN_ROOT/pycache"
export XDG_CACHE_HOME="$RUN_ROOT/cache"
export OPENMM_CPU_THREADS=2
export OMP_NUM_THREADS=2
export PATH="/home/users/lz280/miniconda3/envs/confdiff/bin:$PATH"
# Four small PDB downloads, approximately <10 MiB; capped at 20 MiB per file.
# 82 model inputs expected; all intermediate files retained, budget <1 GiB / 4000 files.
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/prep_pc_experimental.py "$RUN_ROOT/runner.py"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/prep_pc_experimental.slurm "$RUN_ROOT/submission.slurm"
python "$RUN_ROOT/runner.py" "$RUN_ROOT"
