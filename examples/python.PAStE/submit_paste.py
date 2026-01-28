#!/bin/bash
#SBATCH --job-name=sample
#SBATCH -p compsci
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ==== Conda环境初始化 ====
source /home/users/lz280/miniconda3/etc/profile.d/conda.sh
conda init
conda activate AmberTools22

python PASTE_6dv2.py