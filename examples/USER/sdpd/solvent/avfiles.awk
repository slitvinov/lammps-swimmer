#! /usr/bin/awk -f
# Average one column in several files.
# One file:
# x1 y1
# Another file
# x2 y2
# awk -f avfiles.awk -v icol=2 <fst file> <second file> gives
# x2 0.5*(y1+y2)

BEGIN {
    if (length(icol)==0) {
	# a column number to average
	icol=1
    }
}

NF{
    rest[FNR, icol] += $(icol)
    n[FNR]+=1

    # to handle files of the different length
    if (FNR>maxfnr) {
	maxfnr=FNR
    }

    $(icol) = "" 
    time[FNR] = $0
}

END {
    for (q=1; q<maxfnr+1; ++q) {
	print time[q], rest[q, icol]/n[q]
    }
}
