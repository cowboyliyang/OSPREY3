#!/bin/bash
#SBATCH --job-name=ps-log10-all-audit
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_all_audit_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_all_audit_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/audit_${SLURM_JOB_ID}
mkdir "$ROOT"
mkdir "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/audit_all_log10_20260927.py "$ROOT/"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
python3 "$ROOT/audit_all_log10_20260927.py" "$ROOT" | tee "$ROOT/summary.txt"
