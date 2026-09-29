#!/bin/bash
#SBATCH --job-name=ps-log10-incremental
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_incremental_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_incremental_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
BASE=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927
export POSTCHECK_PACKAGE="$BASE/incremental_package_${SLURM_JOB_ID}"
mkdir "$POSTCHECK_PACKAGE"
mkdir "$POSTCHECK_PACKAGE/tmp"
export TMPDIR="$POSTCHECK_PACKAGE/tmp" PYTHONDONTWRITEBYTECODE=1
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/{log10_postcheck_20260927.py,log10_postcheck_20260927.slurm,node_log10_followup_20260927.slurm} "$POSTCHECK_PACKAGE/"
cp "${BASH_SOURCE[0]}" "$POSTCHECK_PACKAGE/submitted.slurm"
cd "$POSTCHECK_PACKAGE"
bash -n log10_postcheck_20260927.slurm node_log10_followup_20260927.slurm
sha256sum log10_postcheck_20260927.py log10_postcheck_20260927.slurm node_log10_followup_20260927.slurm > SHA256SUMS
bash "$POSTCHECK_PACKAGE/log10_postcheck_20260927.slurm"
sbatch --parsable --account=grisman --dependency=afterany:12701384:12701385 --kill-on-invalid-dep=yes \
    --export=ALL,POSTCHECK_PACKAGE="$POSTCHECK_PACKAGE" "$POSTCHECK_PACKAGE/log10_postcheck_20260927.slurm" > "$POSTCHECK_PACKAGE/audit_job.txt"
echo REPLACEMENT_POSTCHECK
cat "$POSTCHECK_PACKAGE/audit_job.txt"
