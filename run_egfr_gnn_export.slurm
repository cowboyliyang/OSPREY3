#!/bin/bash
#SBATCH --job-name=egfr_ccd
#SBATCH --output=bench_logs/egfr_ccd_%j.out
#SBATCH --error=bench_logs/egfr_ccd_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=256G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail

echo "========== EGFR CCD Export =========="
echo "Job: $SLURM_JOB_ID, Node: $(hostname), Start: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

cd /home/users/lz280/IdeaProjects/OSPREY3
mkdir -p bench_logs

export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
export GRADLE_USER_HOME="/tmp/${USER}-gradle-${SLURM_JOB_ID}"
PROJECT_CACHE_DIR="/tmp/osprey-project-cache-${SLURM_JOB_ID}"
mkdir -p "$GRADLE_USER_HOME" "$PROJECT_CACHE_DIR"

OUTPUT_DIR="${OUTPUT_DIR:-gnn_data/egfr_5pos_erlotinib}"
NUM_SAMPLES="${NUM_SAMPLES:-200000}"

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" testClasses 2>&1 | tail -5

./gradlew --no-daemon --project-cache-dir "$PROJECT_CACHE_DIR" test \
    --tests "TestBranchMARKStar.exportGNNDataEGFR5pos" \
    --max-workers=1 \
    -Dorg.gradle.jvmargs="-Xmx230g" \
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
    -Dosprey.gnn.numCPUs=104

echo "Export done: $(date)"
echo ""
echo "Protein data:"
ls -lh "$OUTPUT_DIR/protein/"
echo ""
echo "Complex data:"
ls -lh "$OUTPUT_DIR/complex/"
echo ""
echo "=== DONE ==="
