bdChemo <-
function(x0, xtc, z, xt, bk, Niter = c(1e5,2e5,3e5,1e6,1.5e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,1e7,3e7), Nburn = 2e7, thin = 1e4, sa2l = .5, sa2m = 1, at = 10, bt = 1, al = 25, bl = 6, ql = .025, qu = .975, plot.name = NULL, sample.return = F){

	zs = sort(z,index.return=T)
	z1 = zs$x
	xt1 = xt[zs$ix]
	zt = as.data.frame(table(z1))
	nrep = zt[,2]
	n = length(nrep)
	zu = unique(z1)
	d = as.matrix(dist(zu,diag=T,upper=T))
	zstar=seq(zu[1]-0.1,zu[n]+0.1,length.out=100)

	ns = cumsum(nrep)
	ls = sigma20 = var(bk)
	mubk0 = mean(bk)
	nbk = length(bk)
	
	thin0 = thin
	while ((Niter[4] - Niter[1])/(3*thin0) < 5*n){
		thin0 = thin0/10
	}
	K = length(Niter)
	hpara = matrix(0,Niter[K]/thin0,8)
	philambdas = phimus = matrix(0,Niter[K]/thin0,n)
	update = matrix(NA,Niter[K],5)
		
	alambda0 = 0
	amu0 = 0
	llambda = lmu = (zu[n] - zu[1])/(n-1)
	l2lambda0 = l2mu0 = (llambda)^2
	t2lambda0 = 0.1
	t2mu0 = 0.1
	philambda0 = rnorm(n, alambda0)
	phimu0 = rnorm(n, amu0)
	lambda0 = exp(philambda0)
	mu0 = exp(phimu0)
	km0 = rep(kdlmean(x0,lambda0,mu0),times=nrep)
	kv0 = rep(kdlvar(x0,lambda0,mu0),times=nrep)
	ll0 = normappll(xt1,km0,kv0,sigma20,bk,mubk0)
	gl0 = gfunc(d,l2lambda0)
	gcl0 = chol(gl0)
	igcl0 = solve(gcl0)
	gm0 = gfunc(d,l2mu0)
	gcm0 = chol(gm0)
	igcm0 = solve(gcm0)
	
	siter = 1
	scl = gcl0*sqrt(t2lambda0/200)
	scm = gcm0*sqrt(t2mu0/200)

	for (k in 1:4){
		for (iter in siter:Niter[k]){
			tmp = lambdasample(scl,alambda0,t2lambda0,igcl0,mu0,sigma20,philambda0,lambda0,km0,kv0,ll0,x0,xt1,nrep,bk,mubk0)
			philambda0 = tmp$philambda
			lambda0 = tmp$lambda
			km0 = tmp$km
			kv0 = tmp$kv
			ll0 = tmp$ll
			update[iter,1] = tmp$update

			tmp = musample(scm,amu0,t2mu0,igcm0,lambda0,sigma20,phimu0,mu0,km0,kv0,ll0,x0,xt1,nrep,bk,mubk0)
			phimu0 = tmp$phimu
			mu0 = tmp$mu
			km0 = tmp$km
			kv0 = tmp$kv
			ll0 = tmp$ll
			update[iter,2] = tmp$update
		
			alambda0 = alphasample(igcl0,philambda0,t2lambda0,sa2l)
			amu0 = alphasample(igcm0,phimu0,t2mu0,sa2m)

			if (iter %% thin0 == 0){
				idx = iter/thin0
				hpara[idx, 1:2] = c(alambda0,amu0)
				philambdas[idx, ] = philambda0
				phimus[idx, ] = phimu0
			}
		}
		utmp = mean(update[siter:Niter[k],1])
		if (utmp < 0.15) scl = scl*0.3
		if (utmp > 0.25) scl = scl*1.5
		utmp = mean(update[siter:Niter[k],2])
		if (utmp < 0.15) scm = scm*0.3
		if (utmp > 0.25) scm = scm*1.5
		siter = Niter[k]+1
	}
	covl = cov(philambdas[(Niter[k-1]/thin0+1):(Niter[k]/thin0),])
	scl = chol(covl)*2.4/sqrt(n)
	covm = cov(phimus[(Niter[k-1]/thin0+1):(Niter[k]/thin0),])
	scm = chol(covm)*2.4/sqrt(n)

	malambda = mean(hpara[(Niter[k-1]/thin0+1):(Niter[k]/thin0),1])
	mamu = mean(hpara[(Niter[k-1]/thin0+1):(Niter[k]/thin0),2])
	mphilambda = colMeans(philambdas[(Niter[k-1]/thin0+1):(Niter[k]/thin0),])
	mphimu = colMeans(phimus[(Niter[k-1]/thin0+1):(Niter[k]/thin0),])

	ss = (max(xt1)/100)^2
	sll = slm = 0.1
	alambda0 = rnorm(1,malambda, sqrt(sa2l))
	amu0 = rnorm(1,mamu, sqrt(sa2l))
	philambda0 = mphilambda + crossprod(scl,rnorm(n))
	phimu0 = mphimu + crossprod(scm,rnorm(n))
	sigma20 = rexp(1,1/(ls*100))
	t2lambda0 = rexp(1,10)
	t2mu0 = rexp(1,10)
	l2lambda0 = rnorm(1, llambda^2, 0.5)
	l2mu0 = rnorm(1, lmu^2, 0.5)
	while (l2lambda0 < 0) l2lambda0 = rnorm(1, llambda^2, 0.5)
	while (l2mu0 < 0) l2mu0 = rnorm(1, lmu^2, 0.5)
	gl0 = gfunc(d,l2lambda0)
	gcl0 = chol(gl0)
	igcl0 = solve(gcl0)
	gm0 = gfunc(d,l2mu0)
	gcm0 = chol(gm0)
	igcm0 = solve(gcm0)

	lambda0 = exp(philambda0)
	mu0 = exp(phimu0)
	km0 = rep(kdlmean(x0,lambda0,mu0),times = nrep)
	kv0 = rep(kdlvar(x0,lambda0,mu0),times = nrep)
	ll0 = normappll(xt1,km0,kv0,sigma20,bk,mubk0)
	
	for (k in 5:K){
		for (iter in siter:Niter[k]){

			tmp = lambdasample(scl,alambda0,t2lambda0,igcl0,mu0,sigma20,philambda0,lambda0,km0,kv0,ll0,x0,xt1,nrep,bk,mubk0)
			philambda0 = tmp$philambda
			lambda0 = tmp$lambda
			km0 = tmp$km
			kv0 = tmp$kv
			ll0 = tmp$ll
			update[iter,1] = tmp$update

			tmp = musample(scm,amu0,t2mu0,igcm0,lambda0,sigma20,phimu0,mu0,km0,kv0,ll0,x0,xt1,nrep,bk,mubk0)
			phimu0 = tmp$phimu
			mu0 = tmp$mu
			km0 = tmp$km
			kv0 = tmp$kv
			ll0 = tmp$ll
			update[iter,2] = tmp$update

			tmp = sigma2sample(ss,xt1,km0,kv0,sigma20,ll0,ls,bk,mubk0)
			sigma20 = tmp$sigma2
			ll0 = tmp$ll
			update[iter,3] = tmp$update

			tmp = mubksample(xt1,km0,kv0,sigma20,bk,nbk,ll0,mubk0)
			mubk0 = tmp$mubk
			ll0 = tmp$ll
		
			alambda0 = alphasample(igcl0,philambda0,t2lambda0,sa2l)

			t2lambda0 = tau2sample(igcl0,philambda0,alambda0,at,bt)
				
			amu0 = alphasample(igcm0,phimu0,t2mu0,sa2m)

			t2mu0 = tau2sample(igcm0,phimu0,amu0,at,bt)

			tmp = l2sample(igcl0,l2lambda0,philambda0,alambda0,t2lambda0,d,al,bl,sll)
			l2lambda0 = tmp$l2
			igcl0 = tmp$igc
			update[iter,4]=tmp$update

			tmp = l2sample(igcm0,l2mu0,phimu0,amu0,t2mu0,d,al,bl,slm)
			l2mu0 = tmp$l2
			igcm0 = tmp$igc
			update[iter,5]=tmp$update
		
			if (iter %% thin0 == 0){
				idx = iter/thin0
				hpara[idx,] = c(alambda0,amu0,t2lambda0,t2mu0,sigma20,l2lambda0,l2mu0,mubk0)
				philambdas[idx,] = philambda0
				phimus[idx,] = phimu0
				}
		}
		siter = Niter[k]+1

		utmp = mean(update[(Niter[k-1]+1):Niter[k],1])
		if (utmp < 0.15) scl = scl*0.3
		if (utmp > 0.25) scl = scl*1.5
		utmp = mean(update[(Niter[k-1]+1):Niter[k],2])
		if (utmp < 0.15) scm = scm*0.3
		if (utmp > 0.25) scm = scm*1.5

		utmp = mean(update[(Niter[k-1]+1):Niter[k],3])
		if (utmp < 0.4) ss = ss*0.25
		if (utmp > 0.5) ss = ss*2
		utmp = mean(update[(Niter[k-1]+1):Niter[k],4])
		if (utmp < 0.4) sll = sll*0.25
		if (utmp > 0.5) sll = sll*2
		utmp = mean(update[(Niter[k-1]+1):Niter[k],5])
		if (utmp < 0.4) slm = slm*0.25
		if (utmp > 0.5) slm = slm*2
	}

	est = seq(Nburn/thin0+1, Niter[K]/thin0)
	if (thin != thin0) est = seq(Nburn/thin0+1, Niter[K]/thin0, by = 10)
	philambdastar = phimustar = matrix(0,length(est),length(zstar))
	for (i in 1:length(est)){
		philambdastar[i,] = ratestarsample(zstar,zu,philambdas[est[i],],hpara[est[i],1],hpara[est[i],3],hpara[est[i],6])
		phimustar[i,] = ratestarsample(zstar,zu,phimus[est[i],],hpara[est[i],2],hpara[est[i],4],hpara[est[i],7])
	}

	znew = c(zu,zstar)
	idx = order(znew)
	philambdanew = cbind(philambdas[est,],philambdastar)[,idx]
	phimunew = cbind(phimus[est,],phimustar)[,idx]
	ks = matrix(0,length(est),length(znew))
	lambdanew = exp(philambdanew)
	munew = exp(phimunew)
	for (i in 1:length(est)){
		ks[i,] = kdlmean(x0,lambdanew[i,],munew[i,])
	}
	
	zp = znew[idx]	
	ksmean = apply(ks,2,mean)
	ksbound = apply(ks,2,quantile,c(ql,qu))

	lambdam = colMeans(lambdanew)
	lambdab = apply(lambdanew,2,quantile,c(ql,qu))
	mum = colMeans(munew)
	mub = apply(munew,2,quantile,c(ql,qu))

	if (!is.null(plot.name)){
		pdf(paste(plot.name,".pdf",sep=""))
		plot(zp,lambdam,xlab = "z",ylab = expression(lambda),type = "l",ylim = c(min(lambdab[1,]),max(lambdab[2,])), panel.first=polygon(c(zp,rev(zp)), c(lambdab[1,],rev(lambdab[2,])), col="#ccebc5", border=NA),col="green")
		plot(zp,mum,xlab = "z",ylab = expression(mu),type = "l",ylim = c(min(mub[1,]),max(mub[2,])), panel.first=polygon(c(zp,rev(zp)), c(mub[1,],rev(mub[2,])), col="#fbb4ae", border=NA), col="red")
		plot(zp,ksmean,type = "l",xlab = "z",ylab = "Cell count",ylim = c(min(ksbound[1,],xt),min(max(ksbound[2,]),2*max(xt))), panel.first=polygon(c(zp,rev(zp)), c(ksbound[1,],rev(ksbound[2,])), col="#b3cde3", border=NA), col="blue")
		points(z,xt,pch = 20)
		dev.off()
	}

	GI50 = TGI = LC50 = IC50 = rep(0, length(est))
	xtm = max(mean( xt1[z1==max(z1)] ) - mean(hpara[est,8]),0)
	xtc = max(xtc - mean(hpara[est,8]),0)
	zid1 = which(zp == zu[1])
	for (i in 1:length(est)){
		GI50[i] = 10^(gi50(ks[i,], xtc, x0, zp, zid1))
		TGI[i] = 10^(tgi(ks[i,], x0, zp, zid1))
		LC50[i] = 10^(lc50(ks[i,], x0, zp, zid1))
		IC50[i] = 10^(ic50(ks[i,], xtc, xtm, zp, zid1))
	}

	stat = matrix(0,3,4)
	stat[,1] = c(mean(GI50, na.rm=T), quantile(GI50, c(ql, qu), na.rm=T))
	stat[,2] = c(mean(TGI, na.rm=T), quantile(TGI, c(ql, qu), na.rm=T))
	stat[,3] = c(mean(LC50, na.rm=T), quantile(LC50, c(ql, qu), na.rm=T))
	stat[,4] = c(mean(IC50, na.rm=T), quantile(IC50, c(ql, qu), na.rm=T))

	hparas = apply(hpara[est,], 2, function(x) c(mean(x), quantile(x, c(ql, qu))))

	ifelse (sample.return, return(list(hyperparameters = hparas, znew = zp, lambdas = rbind(lambdam, lambdab), mus = rbind(mum, mub), kmean = rbind(ksmean, ksbound), summary = stat, post.lambda = lambdanew, post.mu = munew)), return(list(hyperparameters = hparas, znew = zp, lambdas = rbind(lambdam, lambdab), mus = rbind(mum, mub), kmean = rbind(ksmean, ksbound), summary = stat)))
}
