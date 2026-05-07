#!/bin/bash
cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p gnn/logs

COMMON="--epochs 200 --batch_size 4096 --lr 1e-3 --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5"
PRETRAINED_CX=gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_checkpoint.pt
PRETRAINED_PR=gnn_data/2RL0_all20_4pos_merged/protein/model/gnn_checkpoint.pt

ARCHS="baseline gat dual_enc multihead transformer"

for CONFSPACE in complex protein; do
    DATA="gnn_data/2RL0_all20_4pos_merged/${CONFSPACE}"
    if [ "$CONFSPACE" = "complex" ]; then
        PRETRAINED=$PRETRAINED_CX
    else
        PRETRAINED=$PRETRAINED_PR
    fi

    for ARCH in $ARCHS; do
        sbatch --job-name="st_${ARCH}_${CONFSPACE:0:2}" \
               --partition=compsci-gpu \
               --gres=gpu:a5000:4 \
               --cpus-per-task=8 \
               --mem=64G \
               --time=8:00:00 \
               --output="gnn/logs/abl_${ARCH}_${CONFSPACE}_%j.out" \
               --error="gnn/logs/abl_${ARCH}_${CONFSPACE}_%j.err" \
               --mail-type=END,FAIL \
               --mail-user=lz280@duke.edu \
               --wrap="source /home/users/lz280/miniconda3/etc/profile.d/conda.sh && conda activate confdiff && python3 gnn/train_subtree_ablation.py --arch ${ARCH} --data_dir ${DATA} ${COMMON}"
    done

    # transfer (full fine-tune)
    sbatch --job-name="st_xfer_${CONFSPACE:0:2}" \
           --partition=compsci-gpu \
           --gres=gpu:a5000:4 \
           --cpus-per-task=8 \
           --mem=64G \
           --time=8:00:00 \
           --output="gnn/logs/abl_transfer_${CONFSPACE}_%j.out" \
           --error="gnn/logs/abl_transfer_${CONFSPACE}_%j.err" \
           --mail-type=END,FAIL \
           --mail-user=lz280@duke.edu \
           --wrap="source /home/users/lz280/miniconda3/etc/profile.d/conda.sh && conda activate confdiff && cd /home/users/lz280/IdeaProjects/OSPREY3 && python3 gnn/train_subtree_ablation.py --arch transfer --data_dir ${DATA} ${COMMON} --pretrained_path ${PRETRAINED}"

    # transfer (frozen backbone)
    sbatch --job-name="st_xfrz_${CONFSPACE:0:2}" \
           --partition=compsci-gpu \
           --gres=gpu:a5000:4 \
           --cpus-per-task=8 \
           --mem=64G \
           --time=8:00:00 \
           --output="gnn/logs/abl_transfer_frozen_${CONFSPACE}_%j.out" \
           --error="gnn/logs/abl_transfer_frozen_${CONFSPACE}_%j.err" \
           --mail-type=END,FAIL \
           --mail-user=lz280@duke.edu \
           --wrap="source /home/users/lz280/miniconda3/etc/profile.d/conda.sh && conda activate confdiff && cd /home/users/lz280/IdeaProjects/OSPREY3 && python3 gnn/train_subtree_ablation.py --arch transfer --data_dir ${DATA} ${COMMON} --pretrained_path ${PRETRAINED} --freeze_backbone"
done
