#!/bin/bash
#SBATCH --job-name=markstar_gnn_eval
#SBATCH --output=markstar_gnn_eval_%j.out
#SBATCH --error=markstar_gnn_eval_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "=============================================="
echo "MARK* CCD vs GNN (eval only, epsilon=0.68, 12 sequences)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"

export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

BASE_DIR="gnn_data/2RL0_markstar_8pos"

# Verify models exist
ls -lh "$BASE_DIR/protein/model/gnn_model.onnx"
ls -lh "$BASE_DIR/complex/model/gnn_model.onnx"

export JAVA_TOOL_OPTIONS="\
-Dosprey.gnn.eval.proteinModelPath=$BASE_DIR/protein/model/gnn_model.onnx \
-Dosprey.gnn.eval.complexModelPath=$BASE_DIR/complex/model/gnn_model.onnx"

echo "Compiling..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" compileJava compileTestJava

echo "Running benchmark..."
./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
  --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.benchmarkMARKStarCCDvsGNN" \
  --max-workers=1 \
  -Dorg.gradle.jvmargs="-Xmx200g"

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
