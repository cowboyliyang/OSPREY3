#!/bin/bash
#SBATCH --job-name=ps_frontier_summary
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${BUILD_ROOT:?}"
: "${ARRAY_JOB:?}"
RESULT_ROOT=/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A${ARRAY_JOB}
mkdir -p "$RESULT_ROOT/summary_tmp"
export TMPDIR="$RESULT_ROOT/summary_tmp"
python3 "$BUILD_ROOT/source/slurm/scripts/summarize_frontier_fullnode_20260917.py" "$RESULT_ROOT"
sacct -j "$ARRAY_JOB" --format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList,ReqTRES%80 -X > "$RESULT_ROOT/slurm_accounting.txt"
