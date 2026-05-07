#!/bin/bash
#SBATCH --job-name=egfr_subtree
#SBATCH --output=bench_logs/egfr_subtree_%j.out
#SBATCH --error=bench_logs/egfr_subtree_%j.err
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

echo "========== EGFR Subtree GNN: protein + complex =========="
echo "Job: $SLURM_JOB_ID, Node: $(hostname), Start: $(date)"

cd /home/users/lz280/IdeaProjects/OSPREY3

eval "$(conda shell.bash hook)"
conda activate confdiff
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/egfr_5pos_erlotinib}"

while [ ! -f "$OUTPUT_DIR/complex/confs.csv" ] || [ ! -f "$OUTPUT_DIR/protein/confs.csv" ]; do
    echo "Waiting for export data..."
    sleep 30
done

echo "GPU 0-3: subtree-protein"
CUDA_VISIBLE_DEVICES=0,1,2,3 python3 gnn/train_subtree.py \
  --data_dir "$OUTPUT_DIR/protein" \
  --epochs 200 --batch_size 8192 --lr 1e-3 \
  --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5 &
PID1=$!

echo "GPU 4-7: subtree-complex"
CUDA_VISIBLE_DEVICES=4,5,6,7 python3 gnn/train_subtree.py \
  --data_dir "$OUTPUT_DIR/complex" \
  --epochs 200 --batch_size 8192 --lr 1e-3 \
  --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5 &
PID2=$!

wait $PID1 && echo "subtree-protein OK" || echo "subtree-protein FAIL"
wait $PID2 && echo "subtree-complex OK" || echo "subtree-complex FAIL"
echo "Subtree training done: $(date)"
echo "=== DONE ==="
