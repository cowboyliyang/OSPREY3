#!/bin/bash
#SBATCH --job-name=lute_sample
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=104
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --output=bench_logs/lute_sample_%j.out
#SBATCH --error=bench_logs/lute_sample_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "LUTE-style sampling + GNN training (4pos × 20AA)"
echo "=============================================="
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

cd /home/users/lz280/IdeaProjects/OSPREY3

OUTPUT_DIR="gnn_data/2RL0_all20_4pos_lute"

# ============================================================
# Step 1: Export data with LUTE-style pairwise A* sampling
# ============================================================
echo ""
echo "========== STEP 1: LUTE-style data export =========="

CP=$(cat bench_logs/test_classpath.txt)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

# numSamples = cap (2M should be enough to cover all pairs × 10)
java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx400g \
  -Dosprey.gnn.numSamples=2000000 \
  -Dosprey.gnn.outputDir=$OUTPUT_DIR \
  -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar export4pos

echo ""
echo "Data export complete: $(date)"
ls -lh "$OUTPUT_DIR"/protein/ "$OUTPUT_DIR"/complex/
echo "Protein confs: $(wc -l < "$OUTPUT_DIR/protein/confs.csv")"
echo "Complex confs: $(wc -l < "$OUTPUT_DIR/complex/confs.csv")"

# ============================================================
# Step 2: Train GNN (both protein and complex)
# ============================================================
echo ""
echo "========== STEP 2: Train GNN =========="

eval "$(conda shell.bash hook)"
conda activate confdiff

nvidia-smi || echo "No GPU"

for SPACE in protein complex; do
    echo ""
    echo "--- Training $SPACE ---"
    python3 gnn/train.py \
      --data_dir "$OUTPUT_DIR/$SPACE" \
      --epochs 500 \
      --aa_embed_dim 32 \
      --pos_embed_dim 16 \
      --node_dim 64 \
      --hidden_dim 128 \
      --num_layers 4 \
      --batch_size 512 \
      --lr 1e-3 \
      --val_frac 0.15 \
      --log_every 50 \
      --seed 42
    echo "$SPACE training done: $(date)"
    ls -lh "$OUTPUT_DIR/$SPACE/model/"
done

echo ""
echo "=============================================="
echo "All done: $(date)"
echo "=============================================="
