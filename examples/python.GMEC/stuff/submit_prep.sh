#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --partition=compsci-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a6000:1
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ==== Conda环境初始化 ====
source /home/users/lz280/miniconda3/etc/profile.d/conda.sh
conda init
conda activate AmberTools22

python ./prepComplex.py --pdbname 6dv2.pdb --cld B --clb H --cutoffWM 6 --cutoffN 4 --heapSize 200000