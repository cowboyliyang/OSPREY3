#!/bin/bash
#SBATCH --job-name=gnn_strat
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --gpus=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:10812816
#SBATCH --output=gnn_stratified_%j.out
#SBATCH --error=gnn_stratified_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

cd /home/users/lz280/IdeaProjects/OSPREY3
PYTHON=/home/users/lz280/miniconda3/envs/confdiff/bin/python
DATA=gnn_data/2RL0_all20_4pos

echo "=== +stratified ablation (protein + complex parallel) ==="
echo "Start: $(date)"

CUDA_VISIBLE_DEVICES=0 $PYTHON gnn/train_ablation.py \
    --data_dir ${DATA} --label protein --only "+stratified" \
    --epochs 500 --batch_size 512 --log_every 100 &

CUDA_VISIBLE_DEVICES=1 $PYTHON gnn/train_ablation.py \
    --data_dir ${DATA} --label complex --only "+stratified" \
    --epochs 500 --batch_size 512 --log_every 100 &

wait
echo "All done: $(date)"
