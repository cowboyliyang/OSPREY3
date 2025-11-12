#!/bin/bash

for k in *-MONTAGE/match*-MONTAGE/kstar*/
do
  cd $k
  sbatch *sh
  cd ../../../
done
