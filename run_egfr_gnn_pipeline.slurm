#!/bin/bash
#SBATCH --job-name=egfr_gnn_pipe
#SBATCH --output=bench_logs/egfr_gnn_pipe_%j.out
#SBATCH --error=bench_logs/egfr_gnn_pipe_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "EGFR Erlotinib 5pos GNN: Export + Train (leaf+subtree, protein+complex)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID, Node: $SLURM_NODELIST, Start: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/egfr_5pos_erlotinib}"
NUM_SAMPLES="${NUM_SAMPLES:-200000}"

# ============================================================
# Step 1: Export (writes $OUTPUT_DIR/protein and $OUTPUT_DIR/complex)
# ============================================================
echo ""
echo "========== STEP 1: Data Export (CCD, 64 CPUs) =========="
echo "Start: $(date)"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" testClasses 2>&1 | tail -5

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
    --tests "TestBranchMARKStar.exportGNNDataEGFR5pos" \
    --max-workers=1 \
    -Dorg.gradle.jvmargs="-Xmx200g" \
    -Dosprey.gnn.egfr.proteinPdb=$OUTPUT_DIR/protein.pdb \
    -Dosprey.gnn.egfr.ligandPdb=$OUTPUT_DIR/ligand.h.pdb \
    -Dosprey.gnn.egfr.prepFile=$OUTPUT_DIR/ligand.prepi \
    -Dosprey.gnn.egfr.frcmodFile=$OUTPUT_DIR/ligand.frcmod \
    -Dosprey.gnn.egfr.tcFile=$OUTPUT_DIR/ligand.tc \
    -Dosprey.gnn.egfr.rotFile=$OUTPUT_DIR/ligand.rot \
    -Dosprey.gnn.egfr.ligandResname=AQ4 \
    -Dosprey.gnn.egfr.hotspotResNums=694,766,768,772,773 \
    -Dosprey.gnn.egfr.outputDir=$OUTPUT_DIR \
    -Dosprey.gnn.numSamples=$NUM_SAMPLES \
    -Dosprey.gnn.numCPUs=64

echo "Export done: $(date)"
ls -lh "$OUTPUT_DIR/protein/" "$OUTPUT_DIR/complex/"

# ============================================================
# Step 2: Train 4 models in parallel — 2 GPUs each on the 8×A5000 node
#   GPU 0,1 -> leaf protein
#   GPU 2,3 -> leaf complex
#   GPU 4,5 -> subtree protein
#   GPU 6,7 -> subtree complex
# ============================================================
echo ""
echo "========== STEP 2: Train 4 GNN models (2 A5000 each) =========="
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate confdiff
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

LEAF_LOG=bench_logs/egfr_leaf_${SLURM_JOB_ID}
SUB_LOG=bench_logs/egfr_subtree_${SLURM_JOB_ID}

# --- leaf protein (GPU 0,1) ---
CUDA_VISIBLE_DEVICES=0,1 python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR/protein" \
  --epochs 500 --aa_embed_dim 32 --pos_embed_dim 16 \
  --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 8192 --lr 1e-3 --val_frac 0.15 --log_every 50 \
  --seed 42 --loss huber --huber_delta 1.0 \
  > ${LEAF_LOG}_protein.log 2>&1 &
PID_LP=$!

# --- leaf complex (GPU 2,3) ---
CUDA_VISIBLE_DEVICES=2,3 python3 gnn/train.py \
  --data_dir "$OUTPUT_DIR/complex" \
  --epochs 500 --aa_embed_dim 32 --pos_embed_dim 16 \
  --node_dim 64 --hidden_dim 128 --num_layers 4 \
  --batch_size 8192 --lr 1e-3 --val_frac 0.15 --log_every 50 \
  --seed 42 --loss huber --huber_delta 1.0 \
  > ${LEAF_LOG}_complex.log 2>&1 &
PID_LC=$!

# --- subtree protein (GPU 4,5) ---
CUDA_VISIBLE_DEVICES=4,5 python3 gnn/train_subtree.py \
  --data_dir "$OUTPUT_DIR/protein" \
  --epochs 200 --batch_size 4096 --lr 1e-3 \
  --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5 \
  > ${SUB_LOG}_protein.log 2>&1 &
PID_SP=$!

# --- subtree complex (GPU 6,7) ---
CUDA_VISIBLE_DEVICES=6,7 python3 gnn/train_subtree.py \
  --data_dir "$OUTPUT_DIR/complex" \
  --epochs 200 --batch_size 4096 --lr 1e-3 \
  --node_dim 32 --hidden_dim 64 --num_layers 3 --log_every 5 \
  > ${SUB_LOG}_complex.log 2>&1 &
PID_SC=$!

echo "PIDs: leaf-prot=$PID_LP, leaf-cplx=$PID_LC, subtree-prot=$PID_SP, subtree-cplx=$PID_SC"
echo "Logs: ${LEAF_LOG}_{protein,complex}.log, ${SUB_LOG}_{protein,complex}.log"

wait $PID_LP && echo "leaf-protein OK"   || echo "leaf-protein FAIL"
wait $PID_LC && echo "leaf-complex OK"   || echo "leaf-complex FAIL"
wait $PID_SP && echo "subtree-protein OK" || echo "subtree-protein FAIL"
wait $PID_SC && echo "subtree-complex OK" || echo "subtree-complex FAIL"

echo ""
echo "All training done: $(date)"
ls -lh "$OUTPUT_DIR/protein/model/" 2>/dev/null || true
ls -lh "$OUTPUT_DIR/complex/model/" 2>/dev/null || true

echo "=============================================="
echo "End: $(date)"
echo "=============================================="
