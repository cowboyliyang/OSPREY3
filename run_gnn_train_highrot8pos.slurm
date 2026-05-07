#!/bin/bash
#SBATCH --job-name=gnn8hr
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:titan_v:4
#SBATCH --mem=700G
#SBATCH --time=48:00:00
#SBATCH --output=bench_logs/gnn8hr_%j.out
#SBATCH --error=bench_logs/gnn8hr_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Full pipeline: export GNN data → train all 4 models in parallel (4 GPUs)
# Usage: sbatch --nodelist=jerry2 run_gnn_train_highrot8pos.slurm

NUM_SAMPLES=${NUM_SAMPLES:-300000}
BASE_DIR="${BASE_DIR:-gnn_data/2RL0_highrot_8pos}"
GNN_EPOCHS=${GNN_EPOCHS:-500}
GNN_EMBED=${GNN_EMBED:-64}
GNN_HIDDEN=${GNN_HIDDEN:-128}
GNN_LAYERS=${GNN_LAYERS:-4}
SUBTREE_EPOCHS=${SUBTREE_EPOCHS:-300}

echo "===== GNN Training Pipeline: 8-pos Highrot (4 GPU) ====="
echo "NumSamples: $NUM_SAMPLES, BaseDir: $BASE_DIR"
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU"

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs "$BASE_DIR"

# ============================================================
# Step 1: Export GNN training data (Java, CPU-intensive CCD)
# ============================================================
if [[ -f "$BASE_DIR/protein/confs.csv" && -f "$BASE_DIR/complex/confs.csv" ]]; then
    echo ""
    echo "STEP 1: Skipped (data already exists)"
    ls -lh "$BASE_DIR/protein/confs.csv" "$BASE_DIR/complex/confs.csv"
else
    echo ""
    echo "=========================================="
    echo "STEP 1: Export GNN Data (protein + complex)"
    echo "=========================================="
    echo "Start: $(date)"

    CP=$(cat bench_logs/test_classpath.txt)
    export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

    java -cp "$CP" \
      --add-modules jdk.incubator.foreign \
      -Xmx600g \
      -Dosprey.gnn.numSamples=$NUM_SAMPLES \
      -Dosprey.gnn.outputDir=$BASE_DIR \
      -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
      edu.duke.cs.osprey.markstar.TestBranchMARKStar exportHighrot8pos

    echo "Data export complete: $(date)"
    echo "Protein data:"; ls -lh "$BASE_DIR/protein/" 2>/dev/null
    echo "Complex data:"; ls -lh "$BASE_DIR/complex/" 2>/dev/null
fi

# ============================================================
# Step 2: Train ALL 4 models in parallel (4 GPUs)
#   GPU 0: protein leaf GNN
#   GPU 1: complex leaf GNN
#   GPU 2: protein subtree GNN
#   GPU 3: complex subtree GNN
# ============================================================
echo ""
echo "=========================================="
echo "STEP 2: Train 4 GNN models in parallel (4 GPUs)"
echo "=========================================="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU visible"

# Protein leaf GNN → GPU 0
echo "[GPU 0] Training protein leaf GNN..."
CUDA_VISIBLE_DEVICES=0 python3 gnn/train.py \
  --data_dir "$BASE_DIR/protein" \
  --epochs "$GNN_EPOCHS" \
  --node_dim "$GNN_EMBED" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 &
PID_P_LEAF=$!

# Complex leaf GNN → GPU 1
echo "[GPU 1] Training complex leaf GNN..."
CUDA_VISIBLE_DEVICES=1 python3 gnn/train.py \
  --data_dir "$BASE_DIR/complex" \
  --epochs "$GNN_EPOCHS" \
  --node_dim "$GNN_EMBED" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 &
PID_C_LEAF=$!

# Protein subtree GNN → GPU 2
echo "[GPU 2] Training protein subtree GNN..."
CUDA_VISIBLE_DEVICES=2 python3 gnn/train_subtree.py \
  --data_dir "$BASE_DIR/protein" \
  --epochs "$SUBTREE_EPOCHS" \
  --node_dim "$GNN_EMBED" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 &
PID_P_SUB=$!

# Complex subtree GNN → GPU 3
echo "[GPU 3] Training complex subtree GNN..."
CUDA_VISIBLE_DEVICES=3 python3 gnn/train_subtree.py \
  --data_dir "$BASE_DIR/complex" \
  --epochs "$SUBTREE_EPOCHS" \
  --node_dim "$GNN_EMBED" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 &
PID_C_SUB=$!

echo "PIDs: protein_leaf=$PID_P_LEAF, complex_leaf=$PID_C_LEAF, protein_subtree=$PID_P_SUB, complex_subtree=$PID_C_SUB"

wait $PID_P_LEAF && echo "[GPU 0] Protein leaf done: $(date)" || echo "[GPU 0] Protein leaf FAILED"
wait $PID_C_LEAF && echo "[GPU 1] Complex leaf done: $(date)" || echo "[GPU 1] Complex leaf FAILED"
wait $PID_P_SUB  && echo "[GPU 2] Protein subtree done: $(date)" || echo "[GPU 2] Protein subtree FAILED"
wait $PID_C_SUB  && echo "[GPU 3] Complex subtree done: $(date)" || echo "[GPU 3] Complex subtree FAILED"

conda deactivate 2>/dev/null || true

echo ""
echo "=========================================="
echo "ALL DONE: $(date)"
echo "=========================================="
echo "Models:"
ls -lh "$BASE_DIR"/*/model/gnn_model.onnx 2>/dev/null || echo "  (leaf models missing)"
ls -lh "$BASE_DIR"/*/model_subtree/subtree_model.onnx 2>/dev/null || echo "  (subtree models missing)"
