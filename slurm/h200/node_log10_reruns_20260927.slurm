#!/bin/bash
#SBATCH --job-name=ps-pro6000-logfix
#SBATCH --account=grisman
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:rtx_pro_6000:4
#SBATCH --time=7-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_pro6000_logfix_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_pro6000_logfix_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${BUILD_ROOT:?Validated output-only patch over the original PRO6000 runtime}"
test -f "$BUILD_ROOT/READY"
ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/node_${SLURM_JOB_ID}
mkdir "$ROOT"
cd "$BUILD_ROOT"
sha256sum -c source.sha256 > "$ROOT/source_checksum.txt"
sha256sum -c launch.sha256 > "$ROOT/launch_checksum.txt"
cp "$BUILD_ROOT/design.slurm" "$ROOT/design.slurm"
cp "${BASH_SOURCE[0]}" "$ROOT/node.slurm"
mkdir "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONPYCACHEPREFIX="$ROOT/pycache"
# Queue independently of the original workers. Refresh observations at launch;
# a separate post-original audit handles cases first observed after this scan.
python3 "$BUILD_ROOT/source/slurm/h200/audit_all_log10_20260927.py" "$ROOT" > "$ROOT/audit_summary.txt"
python3 "$BUILD_ROOT/source/slurm/h200/select_all_log10_reruns_20260927.py" "$ROOT" "$BUILD_ROOT" \
    > "$ROOT/selection.txt"
cat "$ROOT/selection.txt"
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
