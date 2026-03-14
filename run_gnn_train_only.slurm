#!/bin/bash
#SBATCH --job-name=gnn_train
#SBATCH --output=gnn_train_%j.out
#SBATCH --error=gnn_train_%j.err
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "GNN Training Only (data already exported)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU visible"

OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/1CC8_8res}"

python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR" \
  --epochs 300 \
  --embed_dim 32 \
  --hidden_dim 64 \
  --num_layers 3 \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
