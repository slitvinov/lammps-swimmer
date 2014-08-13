#! /bin/bash
# x, y, fx, fy, colorgradient_peratom
/home/litvinov/work/lammps-sph/examples/USER/sph/scripts/lammps2punto.sh data-wall/dump*.*  | awk 'NF{print $3, $4, $6, $7, $9, $10, $2, $NF} !NF' > hm
