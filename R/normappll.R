normappll <-
function(xt,km,kv,sigma2,bk,mubk){
sum(dnorm(xt,km+mubk,sqrt(kv+sigma2),log=TRUE)) + sum(dnorm(bk,mubk,sqrt(sigma2),log=TRUE))
}
