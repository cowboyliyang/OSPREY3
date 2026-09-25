#!/bin/bash
#SBATCH --job-name=ps_1gwc_retry_merge
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
export RESULT_ROOT=/usr/xtmp/lz280/packstar_million_census_20260918/1gwc_20260918
export BUILD_ROOT=/usr/xtmp/lz280/packstar_million_census_20260918/resume_build_12634206
export ARRAY_ID=12637406 SYSTEM=1gwc BLOCK=2048
mkdir -p "$RESULT_ROOT/shard_A$ARRAY_ID"
for task in {0..581}; do
    target="$RESULT_ROOT/shard_A$ARRAY_ID/T$task"
    if [[ ! -e "$target" ]]; then
        ln -s "$RESULT_ROOT/shard_A12634381/T$task" "$target"
    fi
done
exec bash /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/merge_million_census_20260918.slurm
