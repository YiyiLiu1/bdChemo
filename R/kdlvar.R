kdlvar <-
function(x0,lambda,mu){
dif = lambda/mu
v = rep(0,length(dif))
i = which(dif==Inf)
if (length(i) != 0){
v[i] = x0*exp(lambda[i]-mu[i])*(exp(lambda[i]-mu[i])-1)
v[-i] = x0*exp(lambda[-i]-mu[-i])*(exp(lambda[-i]-mu[-i])-1)*((dif[-i]+1)/(dif[-i]-1))
} else{
v = x0*exp(lambda-mu)*(exp(lambda-mu)-1)*((dif+1)/(dif-1))
}
return(v)
}
