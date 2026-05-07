#!/bin/bash
#SBATCH --job-name=s6_debug
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=512G
#SBATCH --time=02:00:00
#SBATCH --output=bench_logs/s6_debug_%j.out
#SBATCH --error=bench_logs/s6_debug_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "===== Strategy6 TIMING DEBUG: 8 problem sequences ====="
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU"

cd /home/users/lz280/IdeaProjects/OSPREY3

GNN_BASE_DIR="gnn_data/2RL0_all20_4pos"
PROTEIN_MODEL="$GNN_BASE_DIR/protein/model/gnn_model.onnx"
COMPLEX_MODEL="$GNN_BASE_DIR/complex/model/gnn_model.onnx"

CP=$(cat bench_logs/test_classpath.txt)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx256g \
  -Dosprey.gnn.benchStrategy=strategy6 \
  -Dosprey.gnn.confSpace=highrot \
  -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
  -Dosprey.gnn.numSeqs=20 \
  -Dosprey.gnn.eval.proteinModelPath=$PROTEIN_MODEL \
  -Dosprey.gnn.eval.complexModelPath=$COMPLEX_MODEL \
  -Dosprey.seqFilter=402,407,422,427,442,447,462,467 \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar

echo "Exit: $?, End: $(date)"
