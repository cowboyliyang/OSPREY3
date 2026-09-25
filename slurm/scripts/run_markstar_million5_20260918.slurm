#!/bin/bash
#SBATCH --job-name=markstar_million5
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%A_%a.err
set -euo pipefail
systems=(4wyu 1gwc 1a0r 4wwi 4kt6)
system=${systems[${SLURM_ARRAY_TASK_ID:?}]}
census=/usr/xtmp/lz280/packstar_million_census_20260918/build_12633933
if [[ $SLURM_ARRAY_TASK_ID -ge 2 ]]; then
    census=/usr/xtmp/lz280/packstar_slow_three_z_20260918/build_12634945
fi
test -f "${CONTROL_BUILD:?}/READY"
out="$CONTROL_BUILD/runs/${system}_J$SLURM_JOB_ID"
mkdir -p "$out/tmp"
export TMPDIR="$out/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp "${BASH_SOURCE[0]}" "$out/runner.slurm"
scontrol show job "$SLURM_JOB_ID" > "$out/slurm.txt"
printf '%s\n' "$census" > "$out/census_root.txt"
sha256sum "$census/input/$system.pdb" "$census/screen/$system.tsv" > "$out/input.sha256"
set +e
/usr/bin/time -o "$out/wall.time" -f 'elapsed=%e maxRssKiB=%M' \
 /home/users/lz280/java/jdk-17.0.2+8/bin/java --add-modules jdk.incubator.foreign \
 --add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED \
 --add-opens java.base/java.lang.invoke=ALL-UNNAMED -XX:-UseSuperWord -Xmx12g \
 -XX:ActiveProcessorCount=4 -Dcontrol.cpus=4 "-Djava.io.tmpdir=$TMPDIR" \
 -cp "$CONTROL_BUILD/classes:$(cat "$census/classpath.txt")" \
 edu.duke.cs.osprey.packstar.MarkStarMillionControl "$census" "$out" "$system" \
 > "$out/run.log" 2>&1
status=$?
set -e
printf '%s\n' "$status" > "$out/exit_code"
exit "$status"
