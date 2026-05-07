#!/bin/bash
#SBATCH --job-name=gnn_bench
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-02
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate confdiff

cd /home/users/lz280/IdeaProjects/OSPREY3
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No nvidia-smi"
python gnn/bench_batch.py
