#!/bin/bash
#SBATCH --job-name=bbk8pos
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --output=bench_logs/bbk8pos_%j.out
#SBATCH --error=bench_logs/bbk8pos_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

# Usage (no Gradle, no lock conflicts — uses pre-built classpath):
#   sbatch --nodelist=jerry7 --cpus-per-task=22 --mem=350G --gres="" run_bbk_gnn_8pos_bench.slurm markstar_ccd
#   sbatch --nodelist=jerry7 --cpus-per-task=22 --mem=350G --gres="" run_bbk_gnn_8pos_bench.slurm branch_ccd
#   sbatch run_bbk_gnn_8pos_bench.slurm markstar_s8
#   sbatch run_bbk_gnn_8pos_bench.slurm branch_s8

METHOD=${1:-markstar_ccd}
NUM_SEQS=${NUM_SEQS:-10}
EPSILON=${EPSILON:-0.683}
MAX_MUT=${MAX_MUT:-8}
GPU_BATCH=${GPU_BATCH:-1000}

# GNN model paths (only needed for S8 methods)
GNN_BASE_DIR="${GNN_BASE_DIR:-gnn_data/2RL0_highrot_8pos}"
PROTEIN_MODEL="${PROTEIN_MODEL:-$GNN_BASE_DIR/protein/model/gnn_model.onnx}"
COMPLEX_MODEL="${COMPLEX_MODEL:-$GNN_BASE_DIR/complex/model/gnn_model.onnx}"
PROTEIN_SUBTREE="${PROTEIN_SUBTREE:-$GNN_BASE_DIR/protein/model_subtree/subtree_model.onnx}"
COMPLEX_SUBTREE="${COMPLEX_SUBTREE:-$GNN_BASE_DIR/complex/model_subtree/subtree_model.onnx}"

echo "===== BBK* + GNN 8-pos Highrot Benchmark ====="
echo "Method: $METHOD"
echo "NumSeqs: $NUM_SEQS, Epsilon: $EPSILON, MaxMut: $MAX_MUT"
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU"

cd /home/users/lz280/IdeaProjects/OSPREY3

# Use pre-built classpath (no Gradle lock conflicts — safe to run in parallel)
CP=$(cat bench_logs/test_classpath.txt)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx$((SLURM_MEM_PER_NODE / 2))m \
  -Dosprey.bbkgnn.method=$METHOD \
  -Dosprey.bbkgnn.numSeqs=$NUM_SEQS \
  -Dosprey.bbkgnn.epsilon=$EPSILON \
  -Dosprey.bbkgnn.numCPUs=$SLURM_CPUS_PER_TASK \
  -Dosprey.bbkgnn.maxMut=$MAX_MUT \
  -Dosprey.gnn.gpuBatchSize=$GPU_BATCH \
  -Dosprey.gnn.eval.proteinModelPath=$PROTEIN_MODEL \
  -Dosprey.gnn.eval.complexModelPath=$COMPLEX_MODEL \
  -Dosprey.gnn.eval.proteinSubtreeModelPath=$PROTEIN_SUBTREE \
  -Dosprey.gnn.eval.complexSubtreeModelPath=$COMPLEX_SUBTREE \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar bbk_gnn_bench

echo "Exit: $?, End: $(date)"
