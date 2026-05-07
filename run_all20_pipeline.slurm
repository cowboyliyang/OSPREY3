#!/bin/bash
#SBATCH --job-name=all20_5pos_gnn
#SBATCH --output=all20_5pos_gnn_%j.out
#SBATCH --error=all20_5pos_gnn_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --nodelist=fennario-02
#SBATCH --gres=gpu:a5000:2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "GNN Pipeline: 5 pos × 20 AA = 3.2M sequences"
echo "  Step 1: Export data (CCD, 104 CPUs)"
echo "  Step 2: Train GNN (2 GPUs)"
echo "  Step 3: CP pfunc evaluation"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID | Node: $SLURM_NODELIST"
echo "Start: $(date)"

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

BASE_DIR="gnn_data/2RL0_all20_5pos"
NUM_SAMPLES=200000

# ============================================================
# Step 1: Export
# ============================================================
echo ""
echo "===== STEP 1: Export (${NUM_SAMPLES} samples) ====="
echo "Start: $(date)"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.numSamples=$NUM_SAMPLES \
-Dosprey.gnn.outputDir=$BASE_DIR \
-Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava 2>&1 | tail -5

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.exportGNNDataAllMutable5pos" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo "Export done: $(date)"
echo "Protein:"
wc -l "$BASE_DIR/protein/"*.csv
echo "Complex:"
wc -l "$BASE_DIR/complex/"*.csv

# ============================================================
# Step 2: Train (2 GPUs in parallel)
# ============================================================
echo ""
echo "===== STEP 2: Train GNN ====="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff
nvidia-smi --query-gpu=name --format=csv,noheader

CUDA_VISIBLE_DEVICES=0 python3 gnn/train.py \
  --data_dir "$BASE_DIR/protein" \
  --epochs 500 --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 512 --lr 1e-3 --val_frac 0.15 --seed 42 --log_every 20 &
PID_P=$!

CUDA_VISIBLE_DEVICES=1 python3 gnn/train.py \
  --data_dir "$BASE_DIR/complex" \
  --epochs 500 --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 512 --lr 1e-3 --val_frac 0.15 --seed 42 --log_every 20 &
PID_C=$!

wait $PID_P && echo "Protein model done: $(date)"
wait $PID_C && echo "Complex model done: $(date)"

# ============================================================
# Step 3: CP evaluation
# ============================================================
echo ""
echo "===== STEP 3: CP Pfunc Evaluation ====="
echo "Start: $(date)"

python3 gnn/cp_pfunc_eval.py --base_dir "$BASE_DIR"

echo ""
echo "=============================================="
echo "All done: $(date)"
echo "=============================================="
