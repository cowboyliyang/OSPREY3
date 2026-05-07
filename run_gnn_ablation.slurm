#!/bin/bash
#SBATCH --job-name=gnn_abl
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --nodelist=fennario-02
#SBATCH --gpus=8
#SBATCH --cpus-per-task=104
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --output=gnn_ablation_%j.out
#SBATCH --error=gnn_ablation_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -e
cd /home/users/lz280/IdeaProjects/OSPREY3

PYTHON=/home/users/lz280/miniconda3/envs/confdiff/bin/python
DATA_200K=gnn_data/2RL0_all20_4pos
DATA_400K=gnn_data/2RL0_all20_4pos_400k

echo "=============================================="
echo "GNN Ablation Study — 8 GPU parallel"
echo "Job: $SLURM_JOB_ID | Node: $SLURM_NODELIST"
echo "Start: $(date)"
echo "=============================================="

# ============================================
# Step 1a: Re-export 200K data (with E1 one-body)
# ============================================
echo ""
echo "===== STEP 1a: Re-export 200K (with E1) ====="

./gradlew compileTestJava 2>&1 | tail -3

JAVA_TOOL_OPTIONS="-Dosprey.gnn.numSamples=200000 -Dosprey.gnn.outputDir=${DATA_200K} -Dosprey.gnn.numCPUs=104" \
  ./gradlew test --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.exportGNNDataAllMutable" 2>&1 | tail -5

echo "200K done: $(date)"
wc -l ${DATA_200K}/protein/onebody_energies.csv ${DATA_200K}/complex/onebody_energies.csv 2>/dev/null || echo "WARNING: onebody not found"

# ============================================
# Step 1b: Export 400K data
# ============================================
echo ""
echo "===== STEP 1b: Export 400K ====="

JAVA_TOOL_OPTIONS="-Dosprey.gnn.numSamples=400000 -Dosprey.gnn.outputDir=${DATA_400K} -Dosprey.gnn.numCPUs=104" \
  ./gradlew test --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.exportGNNDataAllMutable" 2>&1 | tail -5

echo "400K done: $(date)"

# ============================================
# Step 2: Train ALL ablations in parallel (8 GPU × 2 labels = 16 runs)
#   Round 1: 8 variants × protein (8 GPUs parallel)
#   Round 2: 8 variants × complex (8 GPUs parallel)
# ============================================
echo ""
echo "===== STEP 2: Parallel ablation training ====="

VARIANTS="baseline +E1 +pairsum +emat_total +bw_loss +all_feat larger 2x_data"

for LABEL in protein complex; do
    echo ""
    echo "====== ${LABEL}: launching 8 variants on 8 GPUs ======"
    echo "Start: $(date)"

    GPU=0
    PIDS=""
    for V in $VARIANTS; do
        LOG_DIR="${DATA_200K}/ablation/${LABEL}"
        mkdir -p "${LOG_DIR}"

        if [ "$V" = "2x_data" ]; then
            DATA2X_ARG="--data_dir_2x ${DATA_400K}"
        else
            DATA2X_ARG=""
        fi

        echo "  GPU ${GPU}: ${V}"
        CUDA_VISIBLE_DEVICES=$GPU $PYTHON gnn/train_ablation.py \
            --data_dir ${DATA_200K} \
            ${DATA2X_ARG} \
            --label ${LABEL} \
            --only "${V}" \
            --epochs 500 \
            --batch_size 512 \
            --log_every 100 \
            > "${LOG_DIR}/${V}.log" 2>&1 &

        PIDS="$PIDS $!"
        GPU=$((GPU + 1))
    done

    echo "  Waiting for all 8 ${LABEL} jobs (PIDs:${PIDS})..."
    FAIL=0
    for PID in $PIDS; do
        wait $PID || FAIL=$((FAIL + 1))
    done
    echo "${LABEL} done: $(date), failures: ${FAIL}"

    # Print summary from logs
    echo "--- ${LABEL} results ---"
    for V in $VARIANTS; do
        RMSE=$(grep "best=" "${LOG_DIR}/${V}.log" | tail -1 | grep -oP 'best=\K[0-9.]+' || echo "N/A")
        echo "  ${V}: best_val_rmse=${RMSE}"
    done
done

# ============================================
# Step 3: Collect results
# ============================================
echo ""
echo "===== STEP 3: Collect ablation results ====="
for LABEL in protein complex; do
    echo ""
    echo "--- ${LABEL} ---"
    cat "${DATA_200K}/ablation/${LABEL}/ablation_results.json" 2>/dev/null || echo "(no results json)"
    echo ""
    for V in $VARIANTS; do
        echo "=== ${V} ==="
        grep -E "MAE:|RMSE:|Max|Spearman|reduction" "${DATA_200K}/ablation/${LABEL}/${V}.log" 2>/dev/null || true
    done
done

echo ""
echo "=============================================="
echo "All done: $(date)"
echo "=============================================="
