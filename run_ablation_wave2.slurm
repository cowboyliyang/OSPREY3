#!/bin/bash
#SBATCH --job-name=gnn_abl2
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --nodelist=fennario-02
#SBATCH --gpus=8
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --dependency=afterok:10812816
#SBATCH --output=gnn_ablation2_%j.out
#SBATCH --error=gnn_ablation2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

cd /home/users/lz280/IdeaProjects/OSPREY3
PYTHON=/home/users/lz280/miniconda3/envs/confdiff/bin/python
DATA=gnn_data/2RL0_all20_4pos
COMMON="--data_dir ${DATA} --epochs 500 --batch_size 512 --log_every 100"

echo "=============================================="
echo "GNN Ablation Wave 2 — 5 new variants"
echo "  +stratified, +target_norm, +huber, pairwise_mlp, ensemble5"
echo "Job: $SLURM_JOB_ID | Start: $(date)"
echo "=============================================="

# Wave 2 variants (excluding ensemble, which is special)
VARIANTS="+stratified +target_norm +huber pairwise_mlp"

for LABEL in protein complex; do
    echo ""
    echo "====== ${LABEL}: launching variants ======"
    echo "Start: $(date)"

    # 4 simple variants on GPU 0-3
    GPU=0
    PIDS=""
    LOG_DIR="${DATA}/ablation/${LABEL}"
    mkdir -p "${LOG_DIR}"

    for V in $VARIANTS; do
        echo "  GPU ${GPU}: ${V}"
        CUDA_VISIBLE_DEVICES=$GPU $PYTHON gnn/train_ablation.py \
            ${COMMON} --label ${LABEL} --only "${V}" \
            > "${LOG_DIR}/${V}.log" 2>&1 &
        PIDS="$PIDS $!"
        GPU=$((GPU + 1))
    done

    # ensemble5 on GPU 4 (trains 5 models sequentially, takes longer)
    echo "  GPU ${GPU}: ensemble5 (5 seeds)"
    CUDA_VISIBLE_DEVICES=$GPU $PYTHON gnn/train_ablation.py \
        ${COMMON} --label ${LABEL} --only "ensemble5" \
        > "${LOG_DIR}/ensemble5.log" 2>&1 &
    PIDS="$PIDS $!"

    echo "  Waiting for all ${LABEL} jobs..."
    FAIL=0
    for PID in $PIDS; do
        wait $PID || FAIL=$((FAIL + 1))
    done
    echo "${LABEL} done: $(date), failures: ${FAIL}"

    # Print results
    echo "--- ${LABEL} results ---"
    for V in $VARIANTS ensemble5; do
        RMSE=$(grep "best=" "${LOG_DIR}/${V}.log" 2>/dev/null | tail -1 | grep -oP 'best=\K[0-9.]+' || echo "N/A")
        echo "  ${V}: best_val_rmse=${RMSE}"
    done
done

echo ""
echo "=============================================="
echo "All done: $(date)"
echo "=============================================="
