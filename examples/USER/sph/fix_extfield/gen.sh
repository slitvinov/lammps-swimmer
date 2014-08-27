#!/bin/bash

convert -flip img/blue-noise-test.tif tif.txt
awk -f gen.awk tif.txt > table.in
