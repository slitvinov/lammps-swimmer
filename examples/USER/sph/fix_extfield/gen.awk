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
}


NR==1{
    nx=$2
    ny=$3
    printf "variable bn string \"%i %i %i %s\"\n",
	nx, ny, 1, "table.in" > "in.bn"
    close("in.bn")

    Ntotal =  4000
    xmin= 0
    ymin= 0

    xmax = 1.0
    ymax = ymax + (xmax-xmin)/nx*ny
    Vdomain = (xmax - xmin) * (ymax - ymin)
    m      =   Vdomain/Ntotal

    dx = sqrt(Vdomain/Ntotal)
    delta =   1e-1*dx

    printf "region            box block %g %g %g %g %g %g units box\n",
	xmin, xmax, ymin, ymax, -delta, delta > "in.box"
    close("in.box")

    printf  "variable       dx  equal %g\n", dx > "in.lattice"
    printf "lattice sq ${dx} origin 0.5 0.5 0.0\n" >> "in.lattice"
    printf "variable Ntarget equal %i\n", Ntotal >> "in.lattice"
    printf  "variable sph_mass  equal %g\n", m >> "in.lattice"
    close("in.lattice")

    # make a domain bigger for the extrapolations
    xmin -= delta;     ymin -= delta
    xmax += delta;     ymin += delta


    max_g = -1e12
    min_g = 1e12

    FS="[:, ()]*"
    next
}

{

    id++
    gs_R = max(255 - $3, 1)
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

    printf "(gen.awk) m=%g\n", m  > "/dev/stderr"
    printf "(gen.awk) min_g=%i, max_g=%i\n", min_g, max_g > "/dev/stderr"
    printf "(gen.awk) rho_min=%e, rho_max=%e\n", A*min_g, A*max_g > "/dev/stderr"
    printf "(gen.awk) Ncell=%i nx=%i ny=%i\n", Ncell, nx, ny > "/dev/stderr"
}
