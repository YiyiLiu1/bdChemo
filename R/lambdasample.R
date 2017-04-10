lambdasample <-
function(scl,alambda,t2lambda,igcl,mu,sigma2,philambda0,lambda0,km0,kv0,ll0,x0,xt,nrep,bk,mubk){
	n = length(lambda0)
	philambda1 = philambda0+crossprod(scl,rnorm(n))
	lambda1 = exp(philambda1)
	km1 = rep(kdlmean(x0,lambda1,mu),times=nrep)
	kv1 = rep(kdlvar(x0,lambda1,mu),times=nrep)
	ll1 = normappll(xt,km1,kv1,sigma2,bk,mubk)
	prob = ( ll1-ll0 + sum(crossprod(philambda0-alambda,igcl)^2-crossprod(philambda1-alambda,igcl)^2)/(2*t2lambda) )
	if ( !is.nan(prob) & log(runif(1))<prob ){
		return(list(philambda=philambda1,lambda=lambda1,km=km1,kv=kv1,ll=ll1,update=1))
	} else{
		return(list(philambda=philambda0,lambda=lambda0,km=km0,kv=kv0,ll=ll0,update=0))
	}
}
