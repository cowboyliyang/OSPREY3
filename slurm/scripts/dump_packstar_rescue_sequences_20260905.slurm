#!/bin/bash
# Emit the exact global sequence index used by GenericPDBBench/KStar for each
# original design.  This is a read-only conf-space construction job.
#SBATCH --job-name=ps_rescue_seqmap
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

REPO=/home/users/lz280/OSPREY3-fresh-packstar
BUILD_ROOT=${BUILD_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/current_build}
CP_FILE="$BUILD_ROOT/test_classpath.txt"
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
DESIGN_CSV=${DESIGN_CSV:-/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv}
PDB_ROOT=${PDB_ROOT:-/usr/xtmp/lz280/dance_bench/pdbs_prepped}
OUT_ROOT=${OUT_ROOT:-/usr/xtmp/lz280/packstar_rescue_20260905/sequence_map_${SLURM_JOB_ID}}
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

[ -s "$CP_FILE" ] || { echo "missing classpath: $CP_FILE" >&2; exit 2; }
mkdir -p "$OUT_ROOT"
CP=$(<"$CP_FILE")

echo "sequence_map_job=$SLURM_JOB_ID start=$(date --iso-8601=seconds)"
tail -n +3 "$DESIGN_CSV" | while IFS=, read -r design_id pdb description num_seq num_flex mutable flexible rest; do
    [ -n "$design_id" ] || continue
    output="$OUT_ROOT/${design_id}.tsv"
    pdb_path="$PDB_ROOT/$pdb/$pdb.min.reduce.renum.pdb"
    echo "dump design=$design_id pdb=$pdb output=$output"
    "$JAVA" \
        --add-opens java.base/java.util=ALL-UNNAMED \
        --add-opens java.base/java.lang=ALL-UNNAMED \
        -Xmx8g \
        -Dosprey.bench.method=sequence_dump \
        -Dosprey.bench.numCPUs=1 \
        -Dosprey.bench.designId="$design_id" \
        -Dosprey.bench.pdbPath="$pdb_path" \
        -Dosprey.bench.mutable="$mutable" \
        -Dosprey.bench.flexible="$flexible" \
        -Dosprey.bench.outputDir="$OUT_ROOT" \
        -Dosprey.sequenceDump.output="$output" \
        -cp "$CP" "$MAIN" \
        > "$OUT_ROOT/${design_id}.log" 2>&1
done

printf 'sequenceMapJob\t%s\n' "$SLURM_JOB_ID" > "$OUT_ROOT/manifest.tsv"
printf 'designCsv\t%s\n' "$DESIGN_CSV" >> "$OUT_ROOT/manifest.tsv"
printf 'pdbRoot\t%s\n' "$PDB_ROOT" >> "$OUT_ROOT/manifest.tsv"
printf 'buildRoot\t%s\n' "$BUILD_ROOT" >> "$OUT_ROOT/manifest.tsv"
echo "sequence_map_job=$SLURM_JOB_ID end=$(date --iso-8601=seconds) output=$OUT_ROOT"
