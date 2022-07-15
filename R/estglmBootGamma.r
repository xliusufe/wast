estglmBootGamma <- function(data, family = "gaussian", h = NULL, smooth = "sigmoid", weights = "exponential", alpha = NULL, beta = NULL, gamma = NULL, gamma0 = NULL, B = 1000, maxIter = 100, tol = 0.0001) {

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
	z 	= data$U[,-1]
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
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

	if(is.null(alpha)||is.null(beta)||is.null(gamma)){
		dims 	= c(n, p1, p2, p3, maxIter, type)
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

		htheta 	= fitglm$theta
		tbeta 	= fitglm$beta
		halpha 	= tbeta[1:p1]
		hbeta 	= tbeta[-c(1:p1)]
		hdelta 	= (1+z%*%htheta>0)

		zt 		= (z%*%htheta)/h
		if(smooth == 'sigmoid'){
			sh 	= 1.0/(1.0+exp(-zt))
		}
		else if(smooth == 'pnorm'){
			sh 	= pnorm(zt)
		}
		else{
			sh 	= pnorm(zt) + zt*dnorm(zt)
		}
		w1 		= cbind(tx, x*as.numeric(sh))
		xa 		= tx%*%halpha
		xb 		= x%*%hbeta

		halpha 	= tbeta[1:p1]*scal_tx
		hbeta 	= tbeta[-c(1:p1)]*scal_x
		htheta 	= c(1, htheta[-1])
	}
	else{
		halpha	= alpha
		hbeta	= beta
		htheta 	= gamma

		zt 		= (z%*%htheta[-1])/h
		if(smooth == 'sigmoid'){
			sh 	= 1.0/(1.0+exp(-zt))
		}
		else if(smooth == 'pnorm'){
			sh 	= pnorm(zt)
		}
		else{
			sh 	= pnorm(zt) + zt*dnorm(zt)
		}
		w1 		= cbind(tx, x*as.numeric(sh))
		xa 		= tx%*%halpha
		xb 		= x%*%hbeta
		hdelta 	= (1+z%*%htheta[-1]>0)
	}




	hthetaB	= matrix(0, nrow = p3, ncol = B)
	if(is.null(gamma0)){
		gamma0 	= rep(1, p3)
	}
	else{
		gamma0 = gamma0[-1]
	}
	dims 	= c(n, ncol(z), maxIter, type)
	params 	= c(tol, h)
	for(b in 1:B){
		G = switch(weights,
					'exponential'	= rexp(n),
					'norm'			= 1 + rnorm(n),
					'bernoulli'		= 2*rbinom(n,1,prob=0.5)
					)
		if(family=='gaussian'){
			fit = .Call("_EST_LINEAR_WEIGHT_SMOOTH",
					as.numeric(y),
					as.numeric(z),
					as.numeric(gamma0),
					as.numeric(xa),
					as.numeric(xb),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else if(family == 'binomial'){
			fit = .Call("_EST_LOGISTIC_WEIGHT_SMOOTH",
					as.numeric(y),
					as.numeric(z),
					as.numeric(gamma0),
					as.numeric(xa),
					as.numeric(xb),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else if(family == 'poisson'){
			fit = .Call("_EST_POISSON_WEIGHT_SMOOTH",
					as.numeric(y),
					as.numeric(z),
					as.numeric(gamma0),
					as.numeric(xa),
					as.numeric(xb),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
		}
		else{
			stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
		}
		hthetaB[, b]	= fit$coef

	}
	hsigma2 = sqrt(n)*apply(hthetaB, 1, sd)

	return(list(alpha = halpha, beta = hbeta, gamma = htheta, delta = hdelta, std = hsigma2))
}