#!/bin/bash
#SBATCH --job-name=gnn_train
#SBATCH --output=bench_logs/gnn_train_%j.out
#SBATCH --error=bench_logs/gnn_train_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:8
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

export PYTHONUNBUFFERED=1

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU visible"

DATA_DIR="${DATA_DIR:?ERROR: set DATA_DIR, e.g. gnn_data/2RL0_all20_4pos_lute/protein}"
RESUME_CKPT="$DATA_DIR/model/gnn_checkpoint.pt"

echo "Data dir: $DATA_DIR"
RESUME_FLAG=""
if [ -f "$RESUME_CKPT" ]; then
    echo "Resuming from checkpoint: $RESUME_CKPT"
    RESUME_FLAG="--resume $RESUME_CKPT"
fi
echo ""

python3 gnn/train.py \
  --data_dir "$DATA_DIR" \
  --epochs 500 \
  --aa_embed_dim 32 \
  --pos_embed_dim 16 \
  --node_dim 64 \
  --hidden_dim 128 \
  --num_layers 4 \
  --batch_size 4096 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 50 \
  --seed 42 \
  --loss huber \
  --huber_delta 1.0 \
  $RESUME_FLAG

echo "Training done: $(date)"
ls -lh "$DATA_DIR/model/"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
