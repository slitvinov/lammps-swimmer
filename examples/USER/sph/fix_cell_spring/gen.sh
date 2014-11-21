#!/bin/bash

Ntotal=20000
imgfile=img/blue-noise-test.tif

convert ${imgfile} tif.txt
convert -flip ${imgfile} tif_flip.txt
maxima --very-quiet -r "Ntotal: ${Ntotal}$ batchload(\"gen-multi-resolution.mac\")\$"

awk -v Ntotal=${Ntotal} -f gen.awk tif_flip.txt > table.aux
paste -d ' ' table.aux lvl.dat | awk '{print $1, $2, $3, $4, $8}' > table.in

mkdir -p data
awk -v Ntotal=${Ntotal} -f gen-data.awk tif.txt p.dat > data/sph.data
