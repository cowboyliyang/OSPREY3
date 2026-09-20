#!/bin/bash
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/frontier_review_%j.out
#SBATCH --error=/usr/xtmp/lz280/frontier_review_%j.err
set -euo pipefail
export TMPDIR=/usr/xtmp/lz280
python3 /home/users/lz280/OSPREY3-fresh-packstar/slurm/scripts/plan_remaining18_20260920.py
