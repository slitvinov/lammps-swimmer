#! /usr/bin/awk -f

FNR == 1 {
    flag = 0
}


flag {
    # data format
    # ITEM: ATOMS id type x y z vx vy vz c_rho_peratom
    id = $1
    type = $2
    x =  $3; y=$4; z=$5
    vx = $6; vy=$7; vz=$8
    rho = $9
    print rho >> "atom."id
}

/^ITEM: ATOMS id type/ {
    flag = 1
}
