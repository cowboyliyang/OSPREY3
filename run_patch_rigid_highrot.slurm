#!/bin/bash
#SBATCH --job-name=patch_hr
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=bench_logs/patch_hr_%j.out
#SBATCH --error=bench_logs/patch_hr_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

echo "===== Patch E_rigid for highrot 8pos ====="
echo "Node: $(hostname), CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

cd /home/users/lz280/IdeaProjects/OSPREY3

# Compile and update classpath
echo "Compiling..."
./gradlew exportTestClasspath --no-daemon 2>&1 | tail -5
echo "Compile done: $(date)"

CP=$(cat bench_logs/test_classpath.txt)
echo "Classpath length: ${#CP}"
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}"

java -cp "$CP" \
  --add-modules jdk.incubator.foreign \
  -Xmx48g \
  -Dosprey.gnn.numCPUs=$SLURM_CPUS_PER_TASK \
  edu.duke.cs.osprey.markstar.TestBranchMARKStar patchRigidHighrot8pos

echo "Exit: $?, End: $(date)"
