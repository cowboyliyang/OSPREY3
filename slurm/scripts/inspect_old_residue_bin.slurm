#!/bin/bash
#SBATCH --job-name=packstar_old_ccd_inspect
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/old_ccd_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/old_ccd_%j.err
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=grisman
#SBATCH --account=grisman

set -euo pipefail
repo=/home/users/lz280/OSPREY3-fresh-packstar
out=/usr/xtmp/lz280/packstar_migration/old_residue_bin
mkdir -p "$out"
cd "$repo"
git show HEAD:src/main/resources/gpuKernels/cuda/residueCcd.bin \
    > "$out/residueCcd-old.bin"
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
"$cuda_env/bin/cuobjdump" --dump-resource-usage "$out/residueCcd-old.bin"
"$cuda_env/bin/cuobjdump" --dump-elf "$out/residueCcd-old.bin" \
    | sed -n '1,35p'
