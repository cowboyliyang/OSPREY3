#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
~/dlab/henry/projects/TLDR/CCK/osprey3-3.3/bin/osprey3 kstar --complex-confspace ./complex.ccsx --target-confspace ./target.ccsx --design-confspace ./design.ccsx --ensemble-dir ./ensembles --max-simultaneous-mutations 1000000
rm *confdb
echo "completed"