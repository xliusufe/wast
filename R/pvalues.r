gam.init = function(n.initials,q,Z,lb.quantile,ub.quantile,ss=1){
	if(q==1){
		gamma.initials = matrix(1,n.initials,q+1)
		gamma.initials[,1] = -quantile(Z,seq(lb.quantile,ub.quantile,length=n.initials))
	}else{
		n.initials = n.initials/ss
		out = matrix(rnorm(n.initials*q), nrow = n.initials)
		gamma.initials	= out/sqrt(apply(out^2,1,sum))
		Z.gamma.initials = Z %*% t(gamma.initials)
		ll=round(n.initials/ss)
		qtile = sample(seq(lb.quantile,ub.quantile,length=n.initials),n.initials)
		gamma.initials.1 = sapply(1:n.initials,function(x)return(
							-quantile(Z.gamma.initials[,x-floor((x-0.1)/ll)*ll],qtile[x])
						))

		gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,min),gamma.initials.1-0.001,gamma.initials.1)
		gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,max),gamma.initials.1+0.001,gamma.initials.1)
		gamma.initials.aug=do.call("rbind", rep(list(gamma.initials), ss))
		gamma.initials = cbind(gamma.initials.1,gamma.initials.aug)
	}
	return(gamma.initials)
}

EstTn_sst <- function(data, family = "gaussian", K = 2000L, M=2000L) {
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}
	G   	= matrix(rnorm(M*n), nrow = n, ncol = M)

	scal_x 	= 1/sqrt(colSums(x^2))
	x 		= x*matrix(rep(scal_x,each=n), nrow=n, ncol=p2)
	scal_tx	= 1/sqrt(colSums(tx^2))
	tx		= tx*matrix(rep(scal_tx,each=n), nrow=n, ncol=p1)

	fitglm 	= glm(y~tx-1, family = family)
	resuds 	= residuals(fitglm, type = "response")
	weights = fitglm$weights
	dims 	= c(n, p1, p2, p3, K, M)
	fit 	<- .Call("_SST_GLM_bstr",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(resuds),
					as.numeric(weights),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims))
	pvals 	= fit$pvals
	theta 	= fit$theta

	return(pvals)
}

EstTn_slr <- function(data, family = "gaussian", K = 1000L, M=1000L, isApprox = 0) {
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	scal_x 	= 1/sqrt(colSums(x^2))
	x 		= x*matrix(rep(scal_x,each=n), nrow=n, ncol=p2)
	scal_tx	= 1/sqrt(colSums(tx^2))
	tx		= tx*matrix(rep(scal_tx,each=n), nrow=n, ncol=p1)

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}

	if(family=='gaussian'){
		type = 1
	}
	else if(family == 'binomial'){
		type 	= 2
	}
	else if(family == 'poisson'){
		type 	= 3
	}
	else{
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}

	G   	= matrix(rnorm(M*n), nrow = n, ncol = M)
	maxIter = 50
	tol 	= 0.00001
	dims 	= c(n, p1, p2, p3, K, M, type, maxIter, isApprox)

	fit 	<- .Call("_GLM_SLR_C",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(tol))

	pvals	= fit$pvals[1]

	return(pvals)
}

EstTn_ast <- function(data, family = "gaussian", isBeta = 0, shape1 = 0, shape2 = 1, M=1000L) {
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	if(p3>1) isBeta = 0
	maxIter = 50
	tol 	= 0.00001
	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)
	params 	= c(shape1, shape2, tol)

	scal_x 	= 1/sqrt(colSums(x^2))
	x 		= x*matrix(rep(scal_x,each=n), nrow=n, ncol=p2)
	scal_tx	= 1/sqrt(colSums(tx^2))
	tx		= tx*matrix(rep(scal_tx,each=n), nrow=n, ncol=p1)

	fitglm 	= glm(y~tx-1, family = family)
	resids 	= residuals(fitglm, type = "response")
	alphahat= fitglm$coefficients

	mu 		= tx%*%alphahat
	if(family=='gaussian'){
		type = 1
		sig2 = sqrt(sum(resids^2)/(n-p1))
	}
	else if(family == 'binomial'){
		type = 2
		expx = 1/(1+exp(-mu))
	}
	else if(family == 'poisson'){
		type = 3
		expx = exp(mu)
	}
	else{
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}

	yb = matrix(0, n, M)
	for(k in 1:M){
		if(family=='gaussian'){
			yb[,k] 	= rnorm(n, mean = 0, sd = sig2) + mu
		}
		else if(family == 'binomial'){
			yb[,k] 	= runif(n)<expx
		}
		else if(family == 'poisson'){
			yb[,k] 	= rpois(n, expx)
		}
	}
	dims = c(dims, type)
	fitwast <- .Call("_GLM_WAST",
				as.numeric(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.integer(dims),
				as.numeric(params))

	teststat 	= fitwast$Tn0
	teststat_p 	= fitwast$Tns
	pvals   	= mean(teststat_p > teststat)

	return(pvals)
}

pval <- function(data, family = "gaussian", method = "wast", M=1000, K = 2000){
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	if(method=='wast') {
	   pvals  	= EstTn_ast(data, family = family, M=M)
	}
	else if(method=='sst'){
		pvals  	= EstTn_sst(data, family = family, K = K, M=M)
	}
	else if(method=='slrt'){
		pvals  	= EstTn_slr(data, family = family, K = K, M=M)
	}
	else{
		stop("Method must be one of {'wast', 'sst', and 'slrt'} !")
	}
	return(pvals)
}