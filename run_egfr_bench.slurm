#!/bin/bash
#SBATCH --job-name=egfr_bench
#SBATCH --output=bench_logs/egfr_bench_%j.out
#SBATCH --error=bench_logs/egfr_bench_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

cd /home/users/lz280/IdeaProjects/OSPREY3
eval "$(conda shell.bash hook)"
conda activate confdiff

echo "GPU:"; nvidia-smi --query-gpu=name --format=csv,noheader
echo ""

python3 gnn/bench_egfr.py

echo "=== DONE ==="
