sigma2sample <-
function(ss,xt,km,kv,sigma20,ll0,ls,bk,mubk){
	sigma21 = rnorm(1,sigma20,ss)
	if (sigma21<0){
		return(list(sigma2=sigma20,ll=ll0,update=0))
	} else{
		ll1 = normappll(xt,km,kv,sigma21,bk,mubk)
		prob = (ll1-ll0 + log(sigma20) - log(sigma21))
		if ( !is.nan(prob) & log(runif(1))<prob ){
			return(list(sigma2=sigma21,ll=ll1,update=1))
		} else{
			return(list(sigma2=sigma20,ll=ll0,update=0))
		}
	}
}
