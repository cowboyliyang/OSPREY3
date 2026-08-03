#!/bin/bash
# Exact non-time-field comparison of two canonical PACK* merged CSVs.

#SBATCH --job-name=packstar_exactcmp
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=${REPO:-/home/users/lz280/IdeaProjects/OSPREY3}
: "${REFERENCE:?REFERENCE canonical merged CSV is required}"
: "${CANDIDATE:?CANDIDATE canonical merged CSV is required}"
: "${REPORT:?REPORT TSV path under /usr/xtmp/lz280 is required}"
DETAIL=${DETAIL:-${REPORT%.tsv}_detail.tsv}

python3 "$REPO/slurm/wrappers/packstar_compare_merged.py" \
    --reference "$REFERENCE" \
    --candidate "$CANDIDATE" \
    --report "$REPORT" \
    --detail "$DETAIL"
