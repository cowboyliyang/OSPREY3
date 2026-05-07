#!/bin/bash
#SBATCH --job-name=gnn_dbg3
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=bench_logs/gnn_dbg3_%j.out
#SBATCH --error=bench_logs/gnn_dbg3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "Node: $(hostname), Date: $(date)"
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null

cd /home/users/lz280/IdeaProjects/OSPREY3

CP=$(cat bench_logs/test_classpath.txt)

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx32g \
  -Dosprey.gnn.benchStrategy=strategy6 \
  -Dosprey.gnn.numSeqs=1 \
  -Dosprey.gnn.eval.proteinModelPath=gnn_data/2RL0_all20_4pos/protein/model/gnn_model.onnx \
  -Dosprey.gnn.eval.complexModelPath=gnn_data/2RL0_all20_4pos/complex/model/gnn_model.onnx \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar

echo "Exit: $?, Done: $(date)"
