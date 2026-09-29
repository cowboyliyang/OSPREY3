#!/bin/bash
#SBATCH --job-name=ps-log10-complete-plan
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/ps_log10_complete_plan_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/ps_log10_complete_plan_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through Slurm}"
BASE=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/build_12706032
REPO=/home/users/lz280/OSPREY3-fresh-packstar
ROOT=/usr/xtmp/lz280/packstar_pro6000_logfix_20260927/complete_${SLURM_JOB_ID}
mkdir "$ROOT"
mkdir "$ROOT/tmp"
export TMPDIR="$ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
test -f "$BASE/READY"
cd "$BASE"
sha256sum -c source.sha256 > "$ROOT/validated_base_checksums.txt"
sha256sum -c launch.sha256 >> "$ROOT/validated_base_checksums.txt"
# No recompilation: retain the exact production class that passed all 39 tests.
cp -a "$BASE/source" "$ROOT/source"
cp "$BASE/"{test_classpath.txt,git_head.txt,design.slurm,base-source.sha256,base-runtime-source.sha256,regression.txt,output-only.patch} "$ROOT/"
cp "$REPO/slurm/h200/"{audit_all_log10_20260927.py,select_all_log10_reruns_20260927.py,node_log10_reruns_20260927.slurm} "$ROOT/source/slurm/h200/"
cp "${BASH_SOURCE[0]}" "$ROOT/submitted.slurm"
python3 "$ROOT/source/slurm/h200/audit_all_log10_20260927.py" "$ROOT" | tee "$ROOT/audit_summary.txt"
python3 "$ROOT/source/slurm/h200/select_all_log10_reruns_20260927.py" "$ROOT" "$ROOT" | tee "$ROOT/selection.txt"
bash -n "$ROOT/design.slurm" "$ROOT/source/slurm/h200/node_log10_reruns_20260927.slurm"
cmp "$BASE/test_classpath.txt" "$ROOT/test_classpath.txt"
cmp "$BASE/design.slurm" "$ROOT/design.slurm"
cd "$ROOT"
find source -type f -print0 | sort -z | xargs -0 sha256sum > source.sha256
sha256sum "$BASE/classes/edu/duke/cs/osprey/kstar/KStarScore.class" >> source.sha256
sha256sum cases.tsv design.slurm test_classpath.txt > launch.sha256
printf 'validated_patch_build\t%s\nselection_preparation_job\t%s\nselection_refresh\tafter_original_jobs_finish\n' "$BASE" "$SLURM_JOB_ID" > build_manifest.tsv
touch READY
echo "READY=$ROOT"
