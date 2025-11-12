#!/bin/csh -f
#SBATCH --output=all-rounds-log.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
#SBATCH -A grisman --partition=grisman

conda run -n AmberTools23 --live-stream python3 -u main.py

echo "completed"