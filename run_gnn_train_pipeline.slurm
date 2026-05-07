#!/bin/bash
#SBATCH --job-name=gnn_train_pipe
#SBATCH --output=gnn_train_pipe_%j.out
#SBATCH --error=gnn_train_pipe_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=800G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:2
#SBATCH --nodelist=fennario-[01-06]
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "GNN Training Pipeline (Export + Train)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "GPUs: 2"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

# ---- Config (override via env) ----
BASE_DIR="${BASE_DIR:-gnn_data/2RL0_all20_11pos}"
NUM_SAMPLES="${NUM_SAMPLES:-300000}"
EXPORT_TEST="${EXPORT_TEST:-edu.duke.cs.osprey.markstar.TestBranchMARKStar.exportGNNDataAllMutable11pos}"

GNN_EPOCHS="${GNN_EPOCHS:-500}"
GNN_NODE_DIM="${GNN_NODE_DIM:-64}"
GNN_HIDDEN="${GNN_HIDDEN:-128}"
GNN_LAYERS="${GNN_LAYERS:-4}"
GNN_BATCH="${GNN_BATCH:-512}"
GNN_LR="${GNN_LR:-1e-3}"
GNN_LOSS="${GNN_LOSS:-huber}"
HUBER_DELTA="${HUBER_DELTA:-1.0}"

echo "Config:"
echo "  BASE_DIR=$BASE_DIR"
echo "  NUM_SAMPLES=$NUM_SAMPLES"
echo "  EXPORT_TEST=$EXPORT_TEST"
echo "  epochs=$GNN_EPOCHS, node_dim=$GNN_NODE_DIM, hidden_dim=$GNN_HIDDEN, layers=$GNN_LAYERS"
echo "  loss=$GNN_LOSS (delta=$HUBER_DELTA)"

# ============================================================
# Step 1: Export data
# ============================================================
echo ""
echo "===== STEP 1: Export GNN Data ====="
echo "Start: $(date)"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.numSamples=$NUM_SAMPLES \
-Dosprey.gnn.outputDir=$BASE_DIR \
-Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava 2>&1 | tail -5

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "$EXPORT_TEST" \
  --max-workers=1 \
  -DtestMaxHeap=600g

echo "Export complete: $(date)"
echo "Protein:"
ls -lh "$BASE_DIR/protein/" || true
echo "Complex:"
ls -lh "$BASE_DIR/complex/" || true

# ============================================================
# Step 2: Train (2 GPUs in parallel)
# ============================================================
echo ""
echo "===== STEP 2: Train GNN ====="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "No GPU visible"

CUDA_VISIBLE_DEVICES=0 python3 gnn/train.py \
  --data_dir "$BASE_DIR/protein" \
  --epochs "$GNN_EPOCHS" \
  --node_dim "$GNN_NODE_DIM" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size "$GNN_BATCH" \
  --lr "$GNN_LR" \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 \
  --loss "$GNN_LOSS" \
  --huber_delta "$HUBER_DELTA" &
PID_P=$!

CUDA_VISIBLE_DEVICES=1 python3 gnn/train.py \
  --data_dir "$BASE_DIR/complex" \
  --epochs "$GNN_EPOCHS" \
  --node_dim "$GNN_NODE_DIM" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size "$GNN_BATCH" \
  --lr "$GNN_LR" \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42 \
  --loss "$GNN_LOSS" \
  --huber_delta "$HUBER_DELTA" &
PID_C=$!

wait $PID_P && echo "Protein model done: $(date)"
wait $PID_C && echo "Complex model done: $(date)"

conda deactivate 2>/dev/null || true

echo ""
echo "=============================================="
echo "Training pipeline complete: $(date)"
echo "=============================================="
