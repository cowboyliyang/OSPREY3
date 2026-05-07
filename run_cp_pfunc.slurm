#!/bin/bash
#SBATCH --job-name=cp_pfunc
#SBATCH --output=cp_pfunc_%j.out
#SBATCH --error=cp_pfunc_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "CP-Guided Pfunc Computation"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

eval "$(conda shell.bash hook)"
conda activate confdiff

python3 gnn/cp_pfunc_eval.py

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
