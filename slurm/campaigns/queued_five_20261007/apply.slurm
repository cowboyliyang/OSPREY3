#!/bin/bash
#SBATCH --job-name=markstar-queue-five
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_queue_five_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_queue_five_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/queued_five_20261007/prepare.py
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/queued_five_20261007/submit.py
bash /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/report.slurm
