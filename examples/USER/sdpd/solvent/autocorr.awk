#! /usr/bin/awk -f
# get autocorrelation function
# takes mdelta option
NF{
    N++
    array[N]=$1 
}


END  {
    for(i=1; i<=N; i++){
	sum += array[i]
    } 
    mean = sum/N
    for (i=1; i<=N; i++) {
	sum2 += (array[i]-mean)^2
    }
    
    if (!length(mdelta)) mdelta=N
    for (delta=0; delta<mdelta; delta++) {
	sum1=0
	for(i=1; i<=N-delta; i++){
	    sum1 += ((array[i+delta]-mean)*(array[i]-mean))
	}
	# returns biased and un-biased estimate
	print sum1/sum2, sum1*N/(sum2*(N-delta))
    }
}
