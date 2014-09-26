# generate initial conditions

BEGIN {
    FS_BAK=FS
    FS="[:,]"
}

NR==1 && NR==FNR {
    nx=$2
    ny=$3
    xmin= 0
    ymin= 0

    xmax = 1.0
    ymax = ymax + (xmax-xmin)/nx*ny
    Vdomain = (xmax - xmin) * (ymax - ymin)
    m      =   Vdomain/Ntotal

    dx = sqrt(Vdomain/Ntotal)
    delta =   1e-1*dx

    # create intput data file for lammps
    printf "LAMMPS data file via write_data, version 8 Aug 2014, timestep = 0\n"
    printf "\n"
    printf "%i atoms\n", Ntotal
    printf "1 atom types\n"
    printf "%-1.16e %-1.16e xlo xhi\n", xmin, xmax
    printf "%-1.16e %-1.16e ylo yhi\n", ymin, ymax
    printf "%-1.16e %-1.16e zlo zhi\n", -delta, delta
    printf "\n"
    printf "Masses\n"
    printf "\n"
    printf "1 %g\n", m
    printf "\n"
    printf "Atoms # meso\n"
    printf "\n"

    FS=FS_BAK
    nextfile
}

{
    if (n<Ntotal) {
	x =$1;    y=$2
	xc=$3;   yc=$4
	rc=$5
	printf "%i 1 0.0000000000000000e+00 %-1.16e 1.0000000000000000e+00 " \
	"%-1.16e %-1.16e 0.0000000000000000e+00 %-1.16e %-1.16e 0.0000000000000000e+00 0 0 0\n", \
	    ++n, rc, xc, yc, x, y
    }

}

END {
    printf "\n"    
    printf "Velocities\n"
    printf "\n"
    for (n=1; n<=Ntotal; n++) {
	printf "%i 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00\n", n
    }
}


