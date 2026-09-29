#!/bin/bash
#SBATCH --job-name=ps-log10-postcheck-setup
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_postcheck_setup_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_postcheck_setup_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/postcheck_package_${SLURM_JOB_ID}
mkdir "$ROOT"
mkdir "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/{log10_postcheck_20260927.py,log10_postcheck_20260927.slurm,node_log10_followup_20260927.slurm} "$ROOT/"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
cd "$ROOT"
python3 -m py_compile log10_postcheck_20260927.py
bash -n log10_postcheck_20260927.slurm node_log10_followup_20260927.slurm
sha256sum log10_postcheck_20260927.py log10_postcheck_20260927.slurm node_log10_followup_20260927.slurm > SHA256SUMS
sbatch --parsable --account=grisman --dependency=afterany:12701384:12701385 --kill-on-invalid-dep=yes \
    --export=ALL,POSTCHECK_PACKAGE="$ROOT" "$ROOT/log10_postcheck_20260927.slurm" > "$ROOT/audit_job.txt"
printf 'POSTCHECK_PACKAGE=%s\n' "$ROOT"
cat "$ROOT/audit_job.txt"
