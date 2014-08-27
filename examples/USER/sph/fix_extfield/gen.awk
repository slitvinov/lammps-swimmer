function idx2x(i) {
    return xmin + (xmax-xmin)*i/nx
}

function idx2y(j) {
    return ymin + (ymax-ymin)*j/ny
}

function max(x, y) {
    return x>y ? x : y
}


BEGIN {
    FS="[:,]"
    delta = 0.01
    xmin= 0                  - delta
    ymin= 0                  - delta

    xmax = 1.0               + delta
    ymax = 1.0               + delta

    Ntotal =  20408

    m      = 1/Ntotal
    Vdomain = (xmax - xmin) * (ymax - ymin)

    max_g = -1e12
    min_g = 1e12

}


NR==1{
    nx=$2
    ny=$3
    printf "pair_style      sph/bn  %i %i %i %s\n", nx, ny, 1, "table.in" > "in.bn"
    close("in.bn")

    FS="[:, ()]*"
    next
}

{

    id++
    gs_R = max(255 - $3, 0)
    s[id] = gs_R
    xidx[id] = $1
    yidx[id] = $2

    gs_sum += gs_R

    if (gs_R>max_g) {max_g=gs_R}
    if (gs_R<min_g) {min_g=gs_R}
}

END {
    A = m*nx*ny*Ntotal/(Vdomain * gs_sum)
    printf "(gen.awk) A       = %e\n", A > "/dev/stderr"
    printf "(gen.awk) Vdomain = %e\n", Vdomain > "/dev/stderr"
    printf "(gen.awk) gs_sum = %e\n", gs_sum > "/dev/stderr"
    Ncell   = id
    
    for (i=1; i<=Ncell; i++) {
	rho_local = A*s[i]
	print  xidx[i], yidx[i], 0, rho_local
	xp = idx2x(xidx[i])
	yp = idx2y(yidx[i])
	print  xp, yp, 0, rho_local > "table.ref"
    }
    close("table.ref")

    printf "(gen.awk) min_g=%i, max_g=%i\n", min_g, max_g > "/dev/stderr"
    printf "(gen.awk) rho_min=%e, rho_max=%e\n", A*min_g, A*max_g > "/dev/stderr"
    printf "(gen.awk) Ncell=%i nx=%i ny=%i\n", Ncell, nx, ny > "/dev/stderr"
}
