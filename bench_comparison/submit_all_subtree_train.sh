#!/bin/bash
# ===================================================================
# Submit subtree GNN training jobs for all PDBs × {protein,complex}.
# Ligand is excluded — S9 doesn't use a ligand subtree model.
# ===================================================================
set -euo pipefail

GNNDIR=/usr/xtmp/lz280/bench_comparison/gnn_models
TRAIN=/home/users/lz280/IdeaProjects/OSPREY3/gnn/train_subtree.py
LOGDIR=/usr/xtmp/lz280/bench_comparison/logs/train_subtree
SPECS=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv

MIN_CONFS=50

mkdir -p "$LOGDIR"

submitted=0
skipped=0
already=0

while IFS=',' read -r did pdb desc nseq nres mutable flexible nconf trank; do
    [[ "$did" == \#* ]] && continue
    [[ -z "$did" ]] && continue

    for confspace in protein complex; do
        data_dir="$GNNDIR/$pdb/$confspace"
        confs_file="$data_dir/confs.csv"

        # Skip if export not done
        [ -f "$confs_file" ] || continue

        nconfs=$(($(wc -l < "$confs_file") - 1))
        if [ "$nconfs" -lt "$MIN_CONFS" ]; then
            skipped=$((skipped + 1))
            continue
        fi

        # Skip if already trained
        if [ -f "$data_dir/model_subtree/subtree_model.onnx" ]; then
            already=$((already + 1))
            continue
        fi

        # Tier by dataset size (subtree needs more data than leaf)
        if [ "$nconfs" -lt 2000 ]; then
            bs=256; ep=400; nd=8; hd=16; nl=2; dp=0.2; lr=5e-4
            tier="tiny"
        elif [ "$nconfs" -lt 20000 ]; then
            bs=1024; ep=300; nd=16; hd=32; nl=2; dp=0.15; lr=1e-3
            tier="small"
        else
            bs=4096; ep=200; nd=32; hd=64; nl=3; dp=0.1; lr=1e-3
            tier="full"
        fi

        JID=$(sbatch \
            --partition=grisman,compsci-gpu --account=grisman \
            --exclude=jerry1 \
            --gres=gpu:a5000:1 \
            --cpus-per-task=4 \
            --mem=32G \
            --time=02:00:00 \
            --job-name="sub_${pdb}_${confspace:0:2}" \
            --output="$LOGDIR/sub_${pdb}_${confspace}_%j.out" \
            --error="$LOGDIR/sub_${pdb}_${confspace}_%j.err" \
            --mail-user=lz280@duke.edu \
            --mail-type=FAIL \
            --wrap "eval \"\$(/home/users/lz280/miniconda3/bin/conda shell.bash hook)\" && \
                conda activate confdiff && \
                python3 $TRAIN \
                    --data_dir $data_dir \
                    --epochs $ep \
                    --batch_size $bs \
                    --lr $lr \
                    --node_dim $nd \
                    --hidden_dim $hd \
                    --num_layers $nl \
                    --dropout $dp \
                    --log_every 10 && \
                echo \"Done: \$(date)\"" \
            | awk '{print $4}')

        printf "%-6s %-8s %6d confs  %-5s bs=%4d  job %s\n" "$pdb" "$confspace" "$nconfs" "$tier" "$bs" "$JID"
        submitted=$((submitted + 1))
    done
done < "$SPECS"

echo ""
echo "Submitted $submitted subtree training jobs."
echo "Already trained:   $already"
echo "Skipped (<$MIN_CONFS confs): $skipped"
echo "Monitor: squeue -u \$USER --name='sub_%'"
echo "Output:  $GNNDIR/<pdb>/<confspace>/model_subtree/subtree_model.onnx"
