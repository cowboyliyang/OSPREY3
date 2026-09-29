#!/bin/bash
#SBATCH --job-name=ps-pro6000-node
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:rtx_pro_6000:4
#SBATCH --time=7-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_node_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_node_%j.err
set -euo pipefail
: "${BUILD_ROOT:?Validated Duke runtime}"
: "${FIRST_INDEX:?First design index}"
: "${LAST_INDEX:?Last design index}"
test -f "$BUILD_ROOT/READY"
ROOT=/usr/xtmp/lz280/packstar_pro6000_20260925/node_${SLURM_JOB_ID}
mkdir "$ROOT"
cp /home/users/lz280/OSPREY3-fresh-packstar/slurm/h200/run_pro6000_20260925.slurm "$ROOT/design.slurm"
cp "${BASH_SOURCE[0]}" "$ROOT/node.slurm"
bash -n "$ROOT/design.slurm"
printf 'index\tstart\tend\texit_code\tstatus\n' > "$ROOT/status.tsv"
for (( index=FIRST_INDEX; index<=LAST_INDEX; index++ )); do
    export SLURM_ARRAY_TASK_ID="$index"
    start=$(date -Is)
    echo "DESIGN_START index=$index time=$start"
    code=0
    timeout --signal=TERM --kill-after=10s 90m bash "$ROOT/design.slurm" \
        > "$ROOT/design_${index}.out" 2> "$ROOT/design_${index}.err" || code=$?
    case "$code" in
        0) status=COMPLETED ;;
        124|137) status=TIMEOUT ;;
        *) status=FAILED_OR_INCOMPLETE ;;
    esac
    end=$(date -Is)
    printf '%s\t%s\t%s\t%s\t%s\n' "$index" "$start" "$end" "$code" "$status" >> "$ROOT/status.tsv"
    echo "DESIGN_END index=$index status=$status exit=$code time=$end"
done
echo 'ALL_DESIGNS_ATTEMPTED'
