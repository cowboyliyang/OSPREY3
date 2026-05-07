#!/bin/bash
#SBATCH --job-name=s6_batch
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=512G
#SBATCH --time=06:00:00
#SBATCH --output=bench_logs/strat_%j.out
#SBATCH --error=bench_logs/strat_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

STRATEGY=${1:-ccd}
CONF_SPACE=${2:-medium}
GNN_BASE_DIR="${GNN_BASE_DIR:-gnn_data/2RL0_all20_4pos}"
PROTEIN_MODEL="${PROTEIN_MODEL:-$GNN_BASE_DIR/protein/model/gnn_model.onnx}"
COMPLEX_MODEL="${COMPLEX_MODEL:-$GNN_BASE_DIR/complex/model/gnn_model.onnx}"
echo "===== Strategy: $STRATEGY, ConfSpace: $CONF_SPACE ====="
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK, Mem: $SLURM_MEM_PER_NODE MB"
echo "Start: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU"
echo "GNN base dir: $GNN_BASE_DIR"

cd /home/users/lz280/IdeaProjects/OSPREY3
CP=$(cat bench_logs/test_classpath.txt)

# Force JVM to see only SLURM-allocated CPUs (test: thread oversubscription hypothesis)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"
echo "JAVA_TOOL_OPTIONS=$JAVA_TOOL_OPTIONS"
echo "JVM availableProcessors will report: $SLURM_CPUS_PER_TASK"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx256g \
  -Dosprey.gnn.benchStrategy=$STRATEGY \
  -Dosprey.gnn.confSpace=$CONF_SPACE \
  -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
  -Dosprey.gnn.numSeqs=20 \
  -Dosprey.gnn.eval.proteinModelPath=$PROTEIN_MODEL \
  -Dosprey.gnn.eval.complexModelPath=$COMPLEX_MODEL \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar

echo "Exit: $?, End: $(date)"
