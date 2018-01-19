l2sample <-
function(igc0,l20,rate,alpha,tau2,d,al,bl,sl){
    l21 = rnorm(1,l20,sl)
    if (l21<0){
        return(list(l2 = l20, igc = igc0, update = 0))
    } else{
        g1 = gfunc(d,l21)
        gc1 = chol(g1)
        igc1 = t(chol(chol2inv(gc1)))
        prob = sum(crossprod(rate-alpha,igc0)^2 - crossprod(rate-alpha,igc1)^2)/(2*tau2) + sum(log(diag(igc1))) - sum(log(diag(igc0))) + dgamma(l21, al, bl, log=T) - dgamma(l20, al, bl, log=T)
        if (!is.na(prob) & log(runif(1)) < prob){
            return(list(l2 = l21, igc = igc1, update = 1))
        } else{
            return(list(l2 = l20, igc = igc0, update = 0))
        }
    }
}
