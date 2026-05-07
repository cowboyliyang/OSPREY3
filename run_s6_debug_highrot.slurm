#!/bin/bash
#SBATCH --job-name=s6_dbg_hr
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=512G
#SBATCH --time=01:00:00
#SBATCH --output=bench_logs/s6_dbg_hr_%j.out
#SBATCH --error=bench_logs/s6_dbg_hr_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "===== Strategy6 debug, highrot ====="
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "No GPU"

cd /home/users/lz280/IdeaProjects/OSPREY3

CP=$(cat bench_logs/test_classpath.txt)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx256g \
  -Dosprey.gnn.benchStrategy=strategy6 \
  -Dosprey.gnn.confSpace=highrot \
  -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
  -Dosprey.gnn.numSeqs=1 \
  -Dosprey.gnn.eval.proteinModelPath=gnn_data/2RL0_all20_4pos/protein/model/gnn_model.onnx \
  -Dosprey.gnn.eval.complexModelPath=gnn_data/2RL0_all20_4pos/complex/model/gnn_model.onnx \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar

echo "Exit: $?, Done: $(date)"
