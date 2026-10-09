#!/bin/bash
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodelist=compsci-cluster-fitz-35
#SBATCH --job-name=markstar-migrate-dcc4
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:25:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_migrate_dcc4_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_migrate_dcc4_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/dcc_coverage_to_compsci_20261009/migrate.py
bash /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/current_inventory_20261006/report.slurm
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/campaigns/dcc_coverage_to_compsci_20261009/verify.py
