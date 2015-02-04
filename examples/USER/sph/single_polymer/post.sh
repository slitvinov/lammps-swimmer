#!/bin/bash

#cd data
#PYTHONPATH=~/work/Pizza.py/src  python ../dump2ensight.py dump.*.dat
#cd ..

rm            data/atom.*
awk -f id.awk data/dump.*.dat

for id in $(seq 200 3000); do
	awk -f ~/google-svn/awk/correlation/autocorr.awk data/atom.${id} | \
	awk '{print NR-1, $1}' > data/corr.${id}
done

awk -v icol=2 -f ~/google-svn/awk/avfiles.awk data/corr.* > data/av.atom
