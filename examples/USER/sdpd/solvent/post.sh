#!/bin/bash

root=$(pwd)

cd data

rm            atom.*
rm            corr.*
awk -f "${root}"/id.awk dump.*.dat
for id in $(seq 200 3000); do
	awk -f "${root}"/autocorr.awk atom.${id} | \
	awk '{print NR-1, $1}' > corr.${id}
done

awk -v icol=2 -f ${root}/avfiles.awk corr.* > av.atom
