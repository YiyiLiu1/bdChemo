musample <-
function(scm,amu,t2mu,igcm,lambda,sigma2,phimu0,mu0,km0,kv0,ll0,x0,xt,nrep,bk,mubk){
	n = length(mu0)
	phimu1 = phimu0+crossprod(scm,rnorm(n))
	mu1 = exp(phimu1)
	km1 = rep(kdlmean(x0,lambda,mu1),times=nrep)
	kv1 = rep(kdlvar(x0,lambda,mu1),times=nrep)
	ll1 = normappll(xt,km1,kv1,sigma2,bk,mubk)
	prob = ( ll1-ll0 +  sum(crossprod(phimu0-amu,igcm)^2-crossprod(phimu1-amu,igcm)^2)/(2*t2mu) )
	if ( !is.nan(prob) & log(runif(1))<prob ){
		return(list(phimu=phimu1,mu=mu1,km=km1,kv=kv1,ll=ll1,update=1))
	} else{
		return(list(phimu=phimu0,mu=mu0,km=km0,kv=kv0,ll=ll0,update=0))
	}
}
