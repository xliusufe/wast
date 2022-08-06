estglmBoot <- function(data, family = "gaussian", h = NULL, smooth = "sigmoid", weights = "exponential", B = 1000, maxIter = 100, tol = 0.0001) {

	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
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

	dims 	= c(n, p1, p2, p3-1, maxIter, type)
	params 	= c(tol, h)
	if(family=='gaussian'){
		fitglm 	= .Call("_EST_LINEAR",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else if(family == 'binomial'){
		fitglm 	= .Call("_EST_LOGISTICR",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else if(family == 'poisson'){
		fitglm 	= .Call("_EST_POISSON",
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

	htheta 	= c(1, fitglm$theta)
	tbeta 	= fitglm$beta
	halpha 	= tbeta[1:p1]*scal_tx
	hbeta 	= tbeta[-c(1:p1)]*scal_x
	hdelta 	= z%*%htheta>0


	halphaB	= matrix(0, nrow = p1, ncol = B)
	hbetaB	= matrix(0, nrow = p2, ncol = B)
	hthetaB	= matrix(0, nrow = p3-1, ncol = B)

	for(b in 1:B){
		G = switch(weights,
					'exponential'	= rexp(n),
					'norm'			= 1 + rnorm(n),
					'bernoulli'		= 2*rbinom(n,1,prob=0.5)
					)

		if(family=='gaussian'){
			fit = .Call("_EST_LINEAR_Boot",
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
			fit = .Call("_EST_LOGISTICR_Boot",
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
			fit = .Call("_EST_POISSON_Boot",
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
		halphaB[, b]	= tbeta[1:p1]*scal_tx
		hbetaB[, b]		= tbeta[-c(1:p1)]*scal_x
		hthetaB[, b]	= fit$theta
	}
	hsigma2 = sqrt(n)*apply(rbind(halphaB, hbetaB, hthetaB), 1, sd)

	return(list(alpha = halpha, beta = hbeta, gamma = htheta, delta = hdelta, std = hsigma2))
}