#!/bin/bash
#SBATCH --job-name=packstar_pc_candidates
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
export FUNCTIONAL_AUDIT_OUT=/usr/xtmp/lz280/packstar_functional_validation/candidates_$SLURM_JOB_ID
mkdir -p "$FUNCTIONAL_AUDIT_OUT/tmp"
export TMPDIR="$FUNCTIONAL_AUDIT_OUT/tmp"
# Inputs: two existing prepared PDBs, <1 MiB total. Outputs: three small files.
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/audit_functional_candidates.py "$FUNCTIONAL_AUDIT_OUT/audit.py"
python3 "$FUNCTIONAL_AUDIT_OUT/audit.py"
