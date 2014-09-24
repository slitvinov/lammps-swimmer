#!/bin/bash

set -e
set -u
# transform lammps dump files into punto format
awk '(NR>1)&&(FNR==1){printf "\n"; f=0} f{print $3, $4, log($9+0.1); print $3+1.1, $4, log($10+0.1)} /ITEM: ATOMS/{f=1}' $* 

