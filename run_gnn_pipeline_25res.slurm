#!/bin/bash
#SBATCH --job-name=gnn_25res
#SBATCH --output=gnn_25res_%j.out
#SBATCH --error=gnn_25res_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "GNN Pipeline: 1CC8, 25 residues (A26-A50)"
echo "  Step 1: Export data (parallel CCD, 104 CPUs)"
echo "  Step 2: Train GNN (GPU)"
echo "  Step 3: CCD vs GNN benchmark"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

# Isolate gradle state
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

NUM_SAMPLES="${NUM_SAMPLES:-50000}"
OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/1CC8_25res}"
GNN_EPOCHS="${GNN_EPOCHS:-500}"
GNN_EMBED="${GNN_EMBED:-64}"
GNN_HIDDEN="${GNN_HIDDEN:-128}"
GNN_LAYERS="${GNN_LAYERS:-4}"
EVAL_CONFS="${EVAL_CONFS:-200}"

echo "Config:"
echo "  numSamples=$NUM_SAMPLES"
echo "  outputDir=$OUTPUT_DIR"
echo "  epochs=$GNN_EPOCHS"
echo "  embed_dim=$GNN_EMBED, hidden_dim=$GNN_HIDDEN, num_layers=$GNN_LAYERS"
echo "  eval_confs=$EVAL_CONFS"

# ============================================================
# Step 1: Export data (parallel CCD)
# ============================================================
echo ""
echo "=========================================="
echo "STEP 1: Data Export (parallel CCD)"
echo "=========================================="
echo "Start: $(date)"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.numSamples=$NUM_SAMPLES \
-Dosprey.gnn.outputDir=$OUTPUT_DIR"

echo "Compiling..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.exportGNNData1CC8_25res" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo ""
echo "Data export complete: $(date)"
ls -lh "$OUTPUT_DIR"/
echo "Confs: $(wc -l < "$OUTPUT_DIR/confs.csv") lines"

# ============================================================
# Step 2: Train GNN
# ============================================================
echo ""
echo "=========================================="
echo "STEP 2: Train GNN"
echo "=========================================="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU visible (will use CPU)"

python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR" \
  --epochs "$GNN_EPOCHS" \
  --embed_dim "$GNN_EMBED" \
  --hidden_dim "$GNN_HIDDEN" \
  --num_layers "$GNN_LAYERS" \
  --batch_size 512 \
  --lr 1e-3 \
  --val_frac 0.15 \
  --log_every 20 \
  --seed 42

echo ""
echo "Training complete: $(date)"
ls -lh "$OUTPUT_DIR/model/"

# ============================================================
# Step 3: CCD vs GNN Benchmark
# ============================================================
echo ""
echo "=========================================="
echo "STEP 3: CCD vs GNN Benchmark"
echo "=========================================="
echo "Start: $(date)"

conda deactivate 2>/dev/null || true

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.eval.numConfs=$EVAL_CONFS \
-Dosprey.gnn.eval.topK=10 \
-Dosprey.gnn.eval.modelPath=$OUTPUT_DIR/model/gnn_model.onnx"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.minimization.BenchmarkGridDPvsCCD.benchmarkCCDvsGNN1CC8_25res" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
