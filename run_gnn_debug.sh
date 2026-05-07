#!/bin/bash
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1
#SBATCH --nodelist=fennario-01
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --job-name=gnn_debug
#SBATCH --output=gnn_debug_%j.out
#SBATCH --mail-user=lz280@duke.edu
#SBATCH --mail-type=END,FAIL

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo "---"

./gradlew test --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.benchmarkGNNStrategies" \
    -Dosprey.gnn.benchStrategy=strategy6 \
    -Dosprey.gnn.numSeqs=1 \
    --info 2>&1 | tee gnn_debug_output.log

echo "---"
echo "Done: $(date)"
