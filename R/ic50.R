ic50 <-
function(kmeans, mineff, maxeff, z, zid1){
	cc = (mineff + maxeff)/2
	i1 = which(kmeans <= cc)
	if (kmeans[zid1] > cc){
		i2 = min(i1[i1 >= zid1])
	} else{
		i2 = min(i1)
	}
	if (i2>1){
		return( (z[i2]+z[i2-1])/2 )
	} else{
		return(z[i2])
	}
}
