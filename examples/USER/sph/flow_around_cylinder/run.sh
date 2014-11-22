#! /bin/bash

dname=data-wall
mkdir -p ${dname}
~/prefix-mpich/bin/mpirun -np 6  ../../../../src/lmp_linux -in flow.lmp -var dname ${dname}

