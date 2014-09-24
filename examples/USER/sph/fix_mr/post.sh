#!/bin/bash

./join.sh         dump.*.dat > join.dat
./lammps2punto.sh dump.*.dat > punto.dat
