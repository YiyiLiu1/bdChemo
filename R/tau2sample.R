tau2sample <-
function(igc,rate,alpha,at,bt){
	n = length(rate)
	return( 1/rgamma(1,n/2+at,sum(crossprod(rate-alpha,igc)^2)/2+bt) )
}
