#!/bin/bash
# ===================================================================
# Submit GNN training jobs for all 38 PDBs × 3 confspaces
# Skips confspaces with 0 or too few conformations.
# Uses fennario A5000 GPUs on grisman partition.
# ===================================================================
set -euo pipefail

GNNDIR=/usr/xtmp/lz280/bench_comparison/gnn_models
TRAIN=/home/users/lz280/IdeaProjects/OSPREY3/gnn/train.py
LOGDIR=/usr/xtmp/lz280/bench_comparison/logs/train
SPECS=/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv

MIN_CONFS=50       # skip confspaces with fewer conformations
EPOCHS=200
BATCH_SIZE=4096
LR=1e-3

mkdir -p "$LOGDIR"

submitted=0
skipped=0

while IFS=',' read -r did pdb desc nseq nres mutable flexible nconf trank; do
    [[ "$did" == \#* ]] && continue
    [[ -z "$did" ]] && continue

    for confspace in protein ligand complex; do
        data_dir="$GNNDIR/$pdb/$confspace"
        confs_file="$data_dir/confs.csv"

        # Skip if export not done
        if [ ! -f "$confs_file" ]; then
            continue
        fi

        # Skip if too few conformations
        nconfs=$(($(wc -l < "$confs_file") - 1))
        if [ "$nconfs" -lt "$MIN_CONFS" ]; then
            skipped=$((skipped + 1))
            continue
        fi

        # Skip if already trained (ONNX model exists)
        if [ -f "$data_dir/model/gnn_model.onnx" ]; then
            continue
        fi

        # Scale model + hyperparams by dataset size to prevent overfitting
        #   < 2k confs  → tiny  (node=8,  hid=16,  L=2, ~2.8k params)
        #   2k-20k      → small (node=16, hid=32,  L=2, ~9k params)
        #   >= 20k      → full  (node=32, hid=64,  L=3, ~46k params)
        if [ "$nconfs" -lt 2000 ]; then
            bs=256; ep=400; nd=8; hd=16; nl=2; dp=0.2; lr=5e-4
            tier="tiny"
        elif [ "$nconfs" -lt 20000 ]; then
            bs=1024; ep=300; nd=16; hd=32; nl=2; dp=0.15; lr=1e-3
            tier="small"
        else
            bs=$BATCH_SIZE; ep=$EPOCHS; nd=32; hd=64; nl=3; dp=0.1; lr=$LR
            tier="full"
        fi

        JID=$(sbatch \
            --partition=grisman,compsci-gpu --account=grisman \
            --exclude=jerry1 \
            --gres=gpu:a5000:1 \
            --cpus-per-task=4 \
            --mem=32G \
            --time=01:00:00 \
            --job-name="gnn_${pdb}_${confspace:0:2}" \
            --output="$LOGDIR/train_${pdb}_${confspace}_%j.out" \
            --error="$LOGDIR/train_${pdb}_${confspace}_%j.err" \
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
echo "Submitted $submitted training jobs, skipped $skipped (< $MIN_CONFS confs)"
echo "Monitor: squeue -u \$USER --name='gnn_%'"
echo "Output: $GNNDIR/<pdb>/<confspace>/model/gnn_model.onnx"
