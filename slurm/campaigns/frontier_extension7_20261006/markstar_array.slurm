#!/bin/bash
#SBATCH --account=grisman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=192G
#SBATCH --time=14-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/markstar_extension_%A_%a.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/markstar_extension_%A_%a.err
set -euo pipefail
: "${MARKSTAR_EXTENSION_PREP:?}"
: "${MARKSTAR_EXTENSION_GROUP:?}"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
python3 "$MARKSTAR_EXTENSION_PREP/markstar_runner.py"
