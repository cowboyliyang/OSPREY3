#!/bin/bash
# Select component-level PACK* failures that the matching MARK* run measured.
# This job only reads prior artifacts and writes a small candidate manifest.
#SBATCH --job-name=ps_rescue_select
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
PACKSTAR_ROOT=${PACKSTAR_ROOT:-/usr/xtmp/lz280/packstar_adaptive_frequency_severity/decomposition_cost_grisman_multigpu4_38_A12494208}
MARKSTAR_ROOT=${MARKSTAR_ROOT:-/usr/xtmp/lz280/bench_comparison/results}
DESIGN_CSV=${DESIGN_CSV:-/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv}
OUT_ROOT=${OUT_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/selection_${SLURM_JOB_ID}}
OUTPUT=${OUTPUT:-$OUT_ROOT/candidates.tsv}

mkdir -p "$OUT_ROOT"
echo "selector_job=$SLURM_JOB_ID start=$(date --iso-8601=seconds)"
echo "packstar_root=$PACKSTAR_ROOT"
echo "markstar_root=$MARKSTAR_ROOT"
echo "design_csv=$DESIGN_CSV"

python3 "$REPO/slurm/scripts/select_packstar_rescue_cases.py" \
    --packstar-root "$PACKSTAR_ROOT" \
    --markstar-root "$MARKSTAR_ROOT" \
    --design-csv "$DESIGN_CSV" \
    --output "$OUTPUT"

echo "selector_job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
