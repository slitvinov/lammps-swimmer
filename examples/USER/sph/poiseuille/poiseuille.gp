
set xlabel "y"
set ylabel "vx"
plot "<awk 'NF==4{print $2, $4}' data/vx.av | sort -g" u ($1*2e-3):2 t "SPH", \
     "<seq 0 1e-5 1e-3 | tr ',' '.' | ./poiseuille.awk -v t=5" w l t "theory"
