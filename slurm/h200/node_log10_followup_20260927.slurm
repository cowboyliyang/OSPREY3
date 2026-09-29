#!/bin/bash
#SBATCH --job-name=ps-pro6000-logfix-extra
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:rtx_pro_6000:4
#SBATCH --time=7-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_logfix_extra_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_logfix_extra_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${BUILD_ROOT:?Validated output-only patch}"
: "${FOLLOWUP_ROOT:?Frozen supplementary case list}"
test -f "$BUILD_ROOT/READY"
ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/node_${SLURM_JOB_ID}
mkdir "$ROOT"
cd "$BUILD_ROOT"
sha256sum -c source.sha256 > "$ROOT/source_checksum.txt"
sha256sum -c launch.sha256 > "$ROOT/launch_checksum.txt"
cd "$FOLLOWUP_ROOT"
sha256sum -c cases.sha256 > "$ROOT/cases_checksum.txt"
cp "$FOLLOWUP_ROOT/new_cases.tsv" "$ROOT/cases.tsv"
cp "$BUILD_ROOT/design.slurm" "$ROOT/design.slurm"
cp "${BASH_SOURCE[0]}" "$ROOT/node.slurm"
printf 'index\tstart\tend\texit_code\tstatus\n' > "$ROOT/status.tsv"
while IFS=$'\t' read -r index design rest; do
    if [[ "$index" == index ]]; then continue; fi
    [[ "$index" =~ ^[0-9]+$ ]]
    export SLURM_ARRAY_TASK_ID="$index"
    start=$(date -Is)
    echo "DESIGN_START index=$index design=$design time=$start"
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
    echo "DESIGN_END index=$index design=$design status=$status exit=$code time=$end"
done < "$ROOT/cases.tsv"
echo 'ALL_DESIGNS_ATTEMPTED'
