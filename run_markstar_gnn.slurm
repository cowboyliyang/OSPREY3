#!/bin/bash
#SBATCH --job-name=markstar_gnn
#SBATCH --output=markstar_gnn_%j.out
#SBATCH --error=markstar_gnn_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "MARK* CCD vs GNN Comparison Pipeline"
echo "  13 positions (2 mutable + 11 flexible), 12 sequences"
echo "  Step 1: Export protein + complex GNN data"
echo "  Step 2: Train protein + complex GNN (parallel)"
echo "  Step 3: Run MARK*+CCD vs MARK*+GNN"
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

# Isolate gradle state
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

NUM_SAMPLES="${NUM_SAMPLES:-200000}"
BASE_DIR="${BASE_DIR:-gnn_data/2RL0_markstar_13pos}"
GNN_EPOCHS="${GNN_EPOCHS:-500}"
GNN_EMBED="${GNN_EMBED:-64}"
GNN_HIDDEN="${GNN_HIDDEN:-128}"
GNN_LAYERS="${GNN_LAYERS:-4}"

echo "Config:"
echo "  numSamples=$NUM_SAMPLES"
echo "  baseDir=$BASE_DIR"
echo "  epochs=$GNN_EPOCHS"
echo "  embed_dim=$GNN_EMBED, hidden_dim=$GNN_HIDDEN, num_layers=$GNN_LAYERS"

# ============================================================
# Step 1: Export protein + complex GNN data
# ============================================================
echo ""
echo "=========================================="
echo "STEP 1: Export GNN Data (protein + complex)"
echo "=========================================="
echo "Start: $(date)"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.numSamples=$NUM_SAMPLES \
-Dosprey.gnn.outputDir=$BASE_DIR \
-Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK"

echo "Compiling..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.exportGNNDataForMARKStar" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo ""
echo "Data export complete: $(date)"
echo "Protein data:"
ls -lh "$BASE_DIR/protein/"
echo "Complex data:"
ls -lh "$BASE_DIR/complex/"

# ============================================================
# Step 2: Train protein + complex GNN in parallel (2 GPUs)
# ============================================================
echo ""
echo "=========================================="
echo "STEP 2: Train Protein + Complex GNN (parallel, 2 GPUs)"
echo "=========================================="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU visible (will use CPU)"

# Train protein GNN on GPU 0
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
PID_PROTEIN=$!

# Train complex GNN on GPU 1
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
PID_COMPLEX=$!

echo "Protein GNN training PID: $PID_PROTEIN (GPU 0)"
echo "Complex GNN training PID: $PID_COMPLEX (GPU 1)"

wait $PID_PROTEIN
echo "Protein GNN training complete: $(date)"
ls -lh "$BASE_DIR/protein/model/"

wait $PID_COMPLEX
echo "Complex GNN training complete: $(date)"
ls -lh "$BASE_DIR/complex/model/"

conda deactivate 2>/dev/null || true

# ============================================================
# Step 3: MARK* CCD vs GNN Benchmark
# ============================================================
echo ""
echo "=========================================="
echo "STEP 3: MARK* CCD vs GNN Benchmark"
echo "=========================================="
echo "Start: $(date)"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.eval.proteinModelPath=$BASE_DIR/protein/model/gnn_model.onnx \
-Dosprey.gnn.eval.complexModelPath=$BASE_DIR/complex/model/gnn_model.onnx"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.benchmarkMARKStarCCDvsGNN" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
