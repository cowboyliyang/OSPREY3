#!/bin/bash
#SBATCH --job-name=subtree_gnn
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:8
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=fennario-06
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=gnn/logs/train_subtree_%j.out
#SBATCH --error=gnn/logs/train_subtree_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

cd /home/users/lz280/IdeaProjects/OSPREY3

source /home/users/lz280/miniconda3/etc/profile.d/conda.sh
conda activate confdiff

echo "Node: $(hostname)"
echo "GPU:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

CUDA_VISIBLE_DEVICES=0,1,2,3 python3 gnn/train_subtree.py \
    --data_dir gnn_data/2RL0_all20_4pos_merged/complex \
    --epochs 200 \
    --batch_size 4096 \
    --lr 1e-3 \
    --node_dim 32 \
    --hidden_dim 64 \
    --num_layers 3 \
    --log_every 5 &

CUDA_VISIBLE_DEVICES=4,5,6,7 python3 gnn/train_subtree.py \
    --data_dir gnn_data/2RL0_all20_4pos_merged/protein \
    --epochs 200 \
    --batch_size 4096 \
    --lr 1e-3 \
    --node_dim 32 \
    --hidden_dim 64 \
    --num_layers 3 \
    --log_every 5 &

wait
