#!/bin/bash
#SBATCH --job-name=packstar_residue_ccd_compile_debug
#SBATCH --output=/usr/xtmp/lz280/packstar_migration/compile_ccd_%j.out
#SBATCH --error=/usr/xtmp/lz280/packstar_migration/compile_ccd_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=grisman
#SBATCH --account=grisman
#SBATCH --gres=gpu:a5000:1

set -euo pipefail
repo=/home/users/lz280/OSPREY3-fresh-packstar
cuda_env=/home/users/lz280/miniconda3/envs/confdiff
export PATH="$cuda_env/bin:$cuda_env/cuda/bin:$PATH"
out=/usr/xtmp/lz280/packstar_migration/residue_ccd_debug
mkdir -p "$out"
cd "$repo/src/main/resources/gpuKernels/cuda"
for kernel in residueCcd residueCcdBatch; do
    "$cuda_env/bin/nvcc" -fatbin \
        -gencode=arch=compute_70,code=sm_70 \
        -gencode=arch=compute_86,code=sm_86 \
        -maxrregcount=64 -Xptxas=-v -lineinfo \
        -keep -keep-dir "$out" "$kernel.cu" -o "$out/$kernel.bin"
done
"$cuda_env/bin/cuobjdump" --dump-resource-usage "$out/residueCcd.bin"
"$cuda_env/bin/cuobjdump" --dump-resource-usage "$out/residueCcdBatch.bin"
