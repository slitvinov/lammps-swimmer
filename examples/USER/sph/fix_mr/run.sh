#!/bin/bash

#./gen.sh
mpirun.mpich -np 2 ../../../../src/lmp_linux -in sph_nb_lattice.lmp
