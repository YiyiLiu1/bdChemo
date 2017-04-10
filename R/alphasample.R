alphasample <-
function(igc,rate,tau2,sa2){
	cs = colSums(tcrossprod(igc)/tau2)
	s = sum(cs)
	return(rnorm(1,sum(cs*rate)/(s+1/sa2),sqrt(1/(s+1/sa2))))
}
