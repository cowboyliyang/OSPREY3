#!/bin/bash
#SBATCH --job-name=ps_1gwc_extra_retry_link
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
root=/usr/xtmp/lz280/packstar_million_census_20260918/1gwc_20260918
mkdir -p "$root/shard_A12637406"
for task in 458 495; do
    test -f "$root/shard_A12637742/T$task/READY"
    test "$(cat "$root/shard_A12637742/T$task/exit_code")" = 0
    ln -s "$root/shard_A12637742/T$task" "$root/shard_A12637406/T$task"
done
