gi50 <-
function(kmeans, control, initial, z, zid1){
	cc = (control + initial)/2
	if (min(kmeans[-(1:(zid1-1))]) > cc){
		return(NA)
	} else{
		i1 = which(kmeans <= cc)
		if (kmeans[zid1] > cc){
			i2 = min(i1[i1 >= zid1])
		} else{
			i2 = min(i1)
		}
	}
	if (i2>1){
		return( (z[i2]+z[i2-1])/2 )
	} else{
		return(z[i2])
	}
}
