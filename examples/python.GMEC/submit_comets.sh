#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a6000:1
#SBATCH --mem=0
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ==== Conda环境初始化 ====
source /home/users/lz280/miniconda3/etc/profile.d/conda.sh
conda init
conda activate AmberTools22

python comets_trypsin.py
