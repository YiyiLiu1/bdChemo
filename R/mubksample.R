mubksample <-
function(xt,km,kv,sigma2,bk,nbk,ll0,mubk0){
	tsq = 1/(nbk/sigma2 + sum(1/(sigma2+kv)))
	nu = (sum(bk)/sigma2 + sum((xt-km)/(sigma2+kv)))*tsq
	mubk1 = rnorm(1,nu,sqrt(tsq))
	ll1 = normappll(xt,km,kv,sigma2,bk,mubk1)
	return(list(mubk=mubk1, ll=ll1))
}
