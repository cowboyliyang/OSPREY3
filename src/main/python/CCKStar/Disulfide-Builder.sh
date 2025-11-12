#!/bin/csh -f
#SBATCH --output=disulfide-builder.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
#SBATCH -A grisman --partition=grisman

conda run -n AmberTools23 --live-stream python3 -u Build_Disulfide_Bond.py

echo "completed"