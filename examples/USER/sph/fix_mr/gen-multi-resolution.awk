BEGIN {
    FS="[:,( ]*"
}

NR>1{
    if ( (prev_y != $2) && (length(prev_y)>0) ) {
	printf "\n"
    }
    prev_y = $2
    printf "%i ", $5
}

END {
    printf "\n"
}
