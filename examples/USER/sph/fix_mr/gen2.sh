#!/bin/bash

Ntotal=20000
convert img/blue-noise-test.tif tif.txt
convert -flip img/blue-noise-test.tif tif_flip.txt

awk -v Ntotal=${Ntotal} -f gen.awk tif_flip.txt > table.in

maxima --very-quiet -r "Ntotal: ${Ntotal}$ batchload(\"gen-multi-resolution.mac\")\$"
awk -v Ntotal=${Ntotal} -f gen-data.awk tif.txt p.dat > data/sph.data
