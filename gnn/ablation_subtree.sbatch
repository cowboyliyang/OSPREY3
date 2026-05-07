#!/bin/bash
#SBATCH --job-name=subtree_abl
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:8
#SBATCH --cpus-per-task=32
#SBATCH --nodelist=fennario-06
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=gnn/logs/ablation_subtree_%j.out
#SBATCH --error=gnn/logs/ablation_subtree_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

cd /home/users/lz280/IdeaProjects/OSPREY3

source /home/users/lz280/miniconda3/etc/profile.d/conda.sh
conda activate confdiff

echo "Node: $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

DATA=gnn_data/2RL0_all20_4pos_merged/complex
PRETRAINED=gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_checkpoint.pt
COMMON="--data_dir $DATA --epochs 200 --batch_size 4096 --lr 1e-3 --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5"

# 6 architectures on 8 GPUs: baseline+gat share GPU 0-1, rest get 1 each
# But simpler: just assign 1 GPU each (model is tiny), run 6 in parallel + 2 spare

CUDA_VISIBLE_DEVICES=0 python3 gnn/train_subtree_ablation.py --arch baseline  $COMMON &
CUDA_VISIBLE_DEVICES=1 python3 gnn/train_subtree_ablation.py --arch gat       $COMMON &
CUDA_VISIBLE_DEVICES=2 python3 gnn/train_subtree_ablation.py --arch dual_enc  $COMMON &
CUDA_VISIBLE_DEVICES=3 python3 gnn/train_subtree_ablation.py --arch multihead $COMMON &
CUDA_VISIBLE_DEVICES=4 python3 gnn/train_subtree_ablation.py --arch transformer $COMMON &
CUDA_VISIBLE_DEVICES=5 python3 gnn/train_subtree_ablation.py --arch transfer  $COMMON --pretrained_path $PRETRAINED &
CUDA_VISIBLE_DEVICES=6 python3 gnn/train_subtree_ablation.py --arch transfer  $COMMON --pretrained_path $PRETRAINED --freeze_backbone &

wait
echo "All ablations done."
