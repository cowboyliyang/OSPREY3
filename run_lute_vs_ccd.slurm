#!/bin/bash
#SBATCH --job-name=lute_ccd
#SBATCH --output=bench_logs/lute_ccd_%j.out
#SBATCH --error=bench_logs/lute_ccd_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=512G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

set -eo pipefail
export PYTHONUNBUFFERED=1

echo "=============================================="
echo "LUTE vs CCD Accuracy Benchmark (4pos × 20AA)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

cd /home/users/lz280/IdeaProjects/OSPREY3

CP=$(cat bench_logs/test_classpath.txt)
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

NUM_TEST_CONFS="${NUM_TEST_CONFS:-500}"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx400g \
  -Dosprey.lute.numCPUs=$SLURM_CPUS_PER_TASK \
  -Dosprey.lute.numTestConfs=$NUM_TEST_CONFS \
  -Dosprey.lute.maxRMSE=0.1 \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar luteVsCCD

echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
