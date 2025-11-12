#!/bin/bash

FOLDERNAME=$1

for k in "$FOLDERNAME"/*-kstar/kstar-*/
do
  cd $k
  sbatch *sh
  cd ../../..
done