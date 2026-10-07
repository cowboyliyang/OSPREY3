#!/bin/bash
#SBATCH --job-name=markstar-local14
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:08:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_local14_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_local14_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/local_concurrency14_20261007/apply.py
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/refresh.py
python3 /home/users/lz280/notes/markstar_status_20261007.py
python3 /home/users/lz280/notes/markstar_plan_review_20261007.py
python3 /home/users/lz280/notes/render_markstar_plan_table_20261007.py
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/export_plan.py
