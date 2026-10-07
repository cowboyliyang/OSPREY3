#!/bin/bash
#SBATCH --job-name=plan-supplement-apply
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:08:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/plan_supplement_apply_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/plan_supplement_apply_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/plan_supplements_20261007/apply.py
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/refresh.py
python3 /home/users/lz280/notes/markstar_status_20261007.py
python3 /home/users/lz280/notes/markstar_plan_review_20261007.py
python3 /home/users/lz280/notes/render_markstar_plan_table_20261007.py
