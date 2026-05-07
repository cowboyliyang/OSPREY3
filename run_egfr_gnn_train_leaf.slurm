#!/bin/bash
#SBATCH --job-name=egfr_leaf
#SBATCH --output=bench_logs/egfr_leaf_%j.out
#SBATCH --error=bench_logs/egfr_leaf_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "========== EGFR Leaf GNN: protein + complex =========="
echo "Job: $SLURM_JOB_ID, Node: $(hostname), Start: $(date)"

cd /home/users/lz280/IdeaProjects/OSPREY3

eval "$(conda shell.bash hook)"
conda activate confdiff
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/egfr_5pos_erlotinib}"

# Wait for export to finish (data might not exist yet if submitted independently)
while [ ! -f "$OUTPUT_DIR/complex/confs.csv" ] || [ ! -f "$OUTPUT_DIR/protein/confs.csv" ]; do
    echo "Waiting for export data..."
    sleep 30
done

echo "GPU 0-3: leaf-protein"
CUDA_VISIBLE_DEVICES=0,1,2,3 python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR/protein" \
  --epochs 500 --aa_embed_dim 32 --pos_embed_dim 16 \
  --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 16384 --lr 1e-3 --val_frac 0.15 --log_every 50 \
  --seed 42 --loss huber --huber_delta 1.0 &
PID1=$!

echo "GPU 4-7: leaf-complex"
CUDA_VISIBLE_DEVICES=4,5,6,7 python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR/complex" \
  --epochs 500 --aa_embed_dim 32 --pos_embed_dim 16 \
  --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 16384 --lr 1e-3 --val_frac 0.15 --log_every 50 \
  --seed 42 --loss huber --huber_delta 1.0 &
PID2=$!

wait $PID1 && echo "leaf-protein OK" || echo "leaf-protein FAIL"
wait $PID2 && echo "leaf-complex OK" || echo "leaf-complex FAIL"
echo "Leaf training done: $(date)"
echo "=== DONE ==="
