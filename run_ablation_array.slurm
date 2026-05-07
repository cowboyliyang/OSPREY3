#!/bin/bash
#SBATCH --job-name=gnn_abl
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-100
#SBATCH --output=/usr/xtmp/lz280/logs/ablation/task_%a.out
#SBATCH --error=/usr/xtmp/lz280/logs/ablation/task_%a.err

cd /home/users/lz280/IdeaProjects/OSPREY3
PYTHON=/home/users/lz280/miniconda3/envs/confdiff/bin/python
DATA=gnn_data/2RL0_all20_4pos

# Read task params from task list
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gnn/ablation_tasks.txt)
VARIANT=$(echo "$LINE" | cut -f1)
SEED=$(echo "$LINE" | cut -f2)
LABEL=$(echo "$LINE" | cut -f3)

OUT_DIR="${DATA}/ablation_v2/${LABEL}/seed_${SEED}"
mkdir -p "$OUT_DIR"

echo "Task ${SLURM_ARRAY_TASK_ID}: variant=${VARIANT} seed=${SEED} label=${LABEL}"
echo "Node: $(hostname) GPU: ${CUDA_VISIBLE_DEVICES:-auto}"
echo "Start: $(date)"

$PYTHON gnn/train_ablation.py \
    --data_dir ${DATA} \
    --label ${LABEL} \
    --only "${VARIANT}" \
    --out_dir "${OUT_DIR}" \
    --seed ${SEED} \
    --epochs 500 \
    --batch_size 512 \
    --log_every 100

echo "Done: $(date)"
