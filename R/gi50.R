gi50 <-
function(kmeans, control, cinitial, z, zid1){
    cc = (1 + control / cinitial)/2
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
        return( ((cc - kmeans[i2-1]) * z[i2] - (cc - kmeans[i2]) * z[i2-1]) / (kmeans[i2] - kmeans[i2-1]) )
    } else{
        return(NA)
    }
}
