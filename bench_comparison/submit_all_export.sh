#!/bin/bash
# ===================================================================
# Submit CCD export jobs for all 38 PDBs across grisman + compsci
# Maximizes CPU usage: bigger systems get more cores
# ===================================================================
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
LOGDIR=$OUTDIR/logs
SPECS=$OUTDIR/design_specs_prepped.csv
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
GNNDIR=$OUTDIR/gnn_models

JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/java
JARGS="--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.invoke=ALL-UNNAMED -Xmx60g -Xms4g"
MAIN=edu.duke.cs.osprey.markstar.bench.GenericPDBBench

mkdir -p "$LOGDIR" "$GNNDIR"

# Build classpath
CP="build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | grep -v "onnxruntime" | sort -u); do
    CP="$CP:$jar"
done
echo "$CP" > "$LOGDIR/.classpath_export.txt"

# Read all designs, assign partition + cores based on size
i=0
while IFS=',' read -r did pdb desc nseq nres mutable flexible nconf trank; do
    [[ "$did" == \#* ]] && continue  # skip comments
    [[ -z "$did" ]] && continue

    pdbpath="$PDBDIR/$pdb/${pdb}.min.reduce.renum.pdb"
    if [ ! -f "$pdbpath" ]; then
        echo "SKIP $did: PDB not ready"
        continue
    fi

    datadir="$GNNDIR/$did"
    mkdir -p "$datadir"

    # Already exported? (need protein + ligand + complex)
    if [ -f "$datadir/protein/confs.csv" ] && [ -f "$datadir/ligand/confs.csv" ] && [ -f "$datadir/complex/confs.csv" ]; then
        echo "SKIP $did: already exported"
        continue
    fi

    # Assign cores based on num_flex_mut residues
    if [ "$nres" -le 5 ]; then
        cpus=16
    elif [ "$nres" -le 9 ]; then
        cpus=24
    else
        cpus=32
    fi

    # Alternate between grisman and compsci for load balancing
    if (( i % 2 == 0 )); then
        SLURM_PART="--partition=grisman --account=grisman"
        # Exclude jerry1 (broken)
        SLURM_EXCL="--exclude=jerry1,grisman-37"
        mem="64G"
    else
        SLURM_PART="--partition=compsci"
        SLURM_EXCL=""
        mem="60G"
    fi

    JID=$(sbatch $SLURM_PART $SLURM_EXCL \
        --job-name=exp_${did} \
        --cpus-per-task=$cpus \
        --mem=$mem \
        --time=12:00:00 \
        --mail-user=lz280@duke.edu \
        --mail-type=FAIL \
        --output=$LOGDIR/export_${did}_%j.out \
        --error=$LOGDIR/export_${did}_%j.err \
        --wrap "cd /home/users/lz280/IdeaProjects/OSPREY3 && echo \"=== $did ($pdb) on \$(hostname), $cpus CPUs, \$(date) ===\" && $JAVA $JARGS \
            -Dosprey.bench.pdbPath=$pdbpath \
            -Dosprey.bench.mutable='$mutable' \
            -Dosprey.bench.flexible='$flexible' \
            -Dosprey.bench.method=export_gnn \
            -Dosprey.bench.designId=$did \
            -Dosprey.bench.outputDir=$OUTDIR/results \
            -Dosprey.bench.numCPUs=$cpus \
            -Dosprey.gnn.outputDir=$datadir \
            -Dosprey.gnn.numSamples=200000 \
            -cp \"\$(cat $LOGDIR/.classpath_export.txt)\" $MAIN 2>&1 && echo \"Done: \$(date)\"" \
        | awk '{print $4}')

    printf "%-6s %-5s %2d res  %3d CPUs  %-8s  job %s\n" "$did" "$pdb" "$nres" "$cpus" \
        "$([ $((i%2)) -eq 0 ] && echo grisman || echo compsci)" "$JID"
    i=$((i+1))

done < "$SPECS"

echo ""
echo "Submitted $i export jobs. Monitor: squeue -u \$USER"
echo "Output: $GNNDIR/<pdb>/{protein,complex}/confs.csv"
