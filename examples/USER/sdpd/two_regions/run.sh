#!/bin/bash

# http://stackoverflow.com/a/6930607
# make the script fail on error
set -e
set -u

lmp=../../../../src/lmp_linux
mpirun=mpirun.mpich
nproc=2

${mpirun} -np ${nproc} ${lmp} -in in.run
