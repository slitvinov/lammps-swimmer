#!/bin/bash

convert -flop img/quad_ramp.jpg quad_flop.jpg
convert +append img/quad_ramp.jpg quad_flop.jpg  quad_long.jpg
convert -flip quad_long.jpg tif.txt

awk -f gen.awk tif.txt > table.in
