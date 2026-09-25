#!/bin/bash
#SBATCH --job-name=pc-cypa-merge
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail

test -n "${SLURM_JOB_ID:-}"
: "${ARRAY_ROOT:?ARRAY_ROOT must be exported}"
MERGE_ROOT=/usr/xtmp/lz280/packstar_functional_validation/cypa_parallel_final_${SLURM_JOB_ID}
mkdir -p "$MERGE_ROOT/output"

for task in 0 1 2 3; do
    ROOT="$ARRAY_ROOT/task_${task}"
    test -s "$ROOT/output/COMPUTATION_FINISHED"
    test -s "$ROOT/output/results.tsv"
    test -s "$ROOT/output/preflight.tsv"
done

cp "$ARRAY_ROOT/task_0/protocol.json" "$MERGE_ROOT/protocol.json"
cp "$ARRAY_ROOT/task_0/input/sources.tsv" "$MERGE_ROOT/sources.tsv"
head -n 1 "$ARRAY_ROOT/task_0/output/preflight.tsv" > "$MERGE_ROOT/output/preflight.tsv"
head -n 1 "$ARRAY_ROOT/task_0/output/results.tsv" > "$MERGE_ROOT/output/results.tsv"
for task in 0 1 2 3; do
    ROOT="$ARRAY_ROOT/task_${task}"
    tail -n +2 "$ROOT/output/preflight.tsv" >> "$MERGE_ROOT/output/preflight.tsv"
    tail -n +2 "$ROOT/output/results.tsv" >> "$MERGE_ROOT/output/results.tsv"
    cp "$ROOT/output/"*_census.tsv "$MERGE_ROOT/output/"
done

printf 'parallel_array_root\t%s\nvariants\t4\nseeds\t42,43,44\ntarget_epsilon\t0.2\n' \
    "$ARRAY_ROOT" > "$MERGE_ROOT/parallel_manifest.tsv"
printf 'Inspect merged results.\n' > "$MERGE_ROOT/output/COMPUTATION_FINISHED"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/summarize_pc_cypa.py "$MERGE_ROOT/summarize.py"
/home/users/lz280/miniconda3/envs/confdiff/bin/python "$MERGE_ROOT/summarize.py" "$MERGE_ROOT" > "$MERGE_ROOT/summary.stdout"
