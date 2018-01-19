ratestarsample <-
function(zstar,z,rate,alpha,tau2,l2){
znew = c(zstar,z)
nstar = length(zstar)
n = length(z)
d = as.matrix(dist(znew,diag=T,upper=T))
k = (gfunc(d,l2) + diag(1e-10,n+nstar))*tau2
mp = solve(k[(nstar+1):(nstar+n),(nstar+1):(nstar+n)],k[(nstar+1):(nstar+n),1:nstar])
S = k[1:nstar,1:nstar]-crossprod(mp,k[(nstar+1):(nstar+n),1:nstar])
return( alpha + crossprod(mp,(rate-alpha)) + crossprod(chol(S),rnorm(nstar)) )
}
