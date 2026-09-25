#!/bin/bash
#SBATCH --job-name=pc-prep-audit
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
export PC_AUDIT_OUT=/usr/xtmp/lz280/packstar_functional_validation/prep_audit_${SLURM_JOB_ID}
mkdir -p "$PC_AUDIT_OUT/tmp"
export TMPDIR="$PC_AUDIT_OUT/tmp"
export OMP_NUM_THREADS=1
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/audit_pc_experimental_prep.py "$PC_AUDIT_OUT/runner.py"
/home/users/lz280/miniconda3/envs/confdiff/bin/python "$PC_AUDIT_OUT/runner.py" \
    /usr/xtmp/lz280/packstar_functional_validation/experimental_12552431 \
    /usr/xtmp/lz280/packstar_functional_validation/experimental_12552449 \
    /usr/xtmp/lz280/packstar_functional_validation/experimental_12552542
