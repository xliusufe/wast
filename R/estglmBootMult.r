estglmBootMult <- function(data, family = "gaussian", ng = 2, h = NULL, smooth = "sigmoid", weights = "exponential", B = 1000, maxIter = 100, tol = 0.0001) {

	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	if(ng<1){
		warning("The number of change-planes must be equal or greater than one.  The default ng = 2 is used!")
		ng = 2
	}
	else if(ng==1){
		warning("The number of multiple change-planes is equal or greater than two. If ng = 1, the single change-plane model is used, which is same as function 'estglmBoot'!")
		fit = estglmBoot(data, family = family, h = h, smooth = smooth, maxIter = maxIter, tol = tol)
		return(fit)
	}
	if(!(weights %in% c('exponential', 'norm','bernoulli'))){
		weights = 'exponential'
		warning("Input weights is not one of {'exponential', 'norm', 'bernoulli'}, The default weight 'exponential' is used!")
	}
	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		smooth = 'sigmoid'
		warning("Input smooth function is not one of {'sigmoid', 'pnorm', 'mixnorm'},  The default smooth function 'sigmoid' is used!")
	}
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))
	h 	= ifelse(is.null(h), sqrt(n)/log(n), 1/h)

	scal_tx	= 1/sqrt(colSums(tx^2))
	tx 		= tx*matrix(rep(scal_tx,each=n), nrow=n, ncol=p1)
	scal_x 	= 1/sqrt(colSums(x^2))
	x 		= x*matrix(rep(scal_x,each=n), nrow=n, ncol=p2)

	type = switch(smooth,
					'sigmoid'	= 1,
					'pnorm'		= 2,
					'mixnorm'	= 3
				)

	abeta0 	= c(rep(0, p1+ng*p2), rep(1,p3-1), seq(-0.5,0.5,length.out = ng))
	dims 	= c(n, p1, p2, p3-1, maxIter, type, ng)
	params 	= c(tol, h)

	if(family=='gaussian'){
		fitglm 	= .Call("_EST_LINEAR_MULTI",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else if(family == 'binomial'){
		fitglm 	= .Call("_EST_LOGISTICR_MULTI",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else if(family == 'poisson'){
		fitglm 	= .Call("_EST_POISSON_MULTI",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else{
		stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	htheta 	= fitglm$theta
	tbeta 	= fitglm$beta
	halpha 	= tbeta[1:p1]*scal_tx
	ha 		= fitglm$ha
	hbeta 	= rep(0, p2*ng)
	hdelta 	= matrix(0, nrow = n,  ncol = ng)
	for(k in 1:ng){
		hbeta[((k-1)*p2+1):(k*p2)] = tbeta[(p1+(k-1)*p2+1):(p1+k*p2)]*scal_x
		hdelta[,k] 	= (z[,1]+z[,-1]%*%htheta>ha[k])
	}
	htheta 	= c(1, htheta)


	halphaB	= matrix(0, nrow = p1, ncol = B)
	hbetaB	= matrix(0, nrow = ng*p2, ncol = B)
	hthetaB	= matrix(0, nrow = p3-1, ncol = B)
	haB		= matrix(0, nrow = ng, ncol = B)

	for(b in 1:B){
		G = switch(weights,
					'exponential'	= rexp(n),
					'norm'			= 1 + rnorm(n),
					'bernoulli'		= 2*rbinom(n,1,prob=0.5)
					)
		if(family=='gaussian'){
			fit = .Call("_EST_LINEAR_MULTI_BOOT",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else if(family == 'binomial'){
			fit = .Call("_EST_LOGISTICR_MULTI_BOOT",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else if(family == 'poisson'){
			fit = .Call("_EST_POISSON_MULTI_BOOT",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else{
			stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
		}
		tbeta 	= fit$beta
		halphaB[, b]= tbeta[1:p1]*scal_tx
		for(k in 1:ng){
			hbetaB[((k-1)*p2+1):(k*p2),b] = tbeta[(p1+(k-1)*p2+1):(p1+k*p2)]*scal_x
		}
		hthetaB[, b]	= fit$theta
		haB[, b]		= fit$ha
	}
	hsigma2 = sqrt(n)*apply(rbind(halphaB, hbetaB, hthetaB, haB), 1, sd)

	return(list(alpha = halpha, beta = hbeta, gamma = htheta, delta = hdelta, ha = ha, std = hsigma2))
}