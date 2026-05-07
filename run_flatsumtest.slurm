#!/bin/bash
#SBATCH --job-name=flatsum_vs_markstar
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=flatsum_bench_%j.out
#SBATCH --error=flatsum_bench_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lz280@duke.edu

cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=== Running testFlatSumVsOriginalMARKStar ==="
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestBranchMARKStar.testFlatSumVsOriginalMARKStar" \
    -DtestMaxHeap=32g \
    --info 2>&1 | tail -200

echo "=== Done ==="
