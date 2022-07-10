estglm <- function(data, family = "gaussian", h = NULL, smooth = "sigmoid", maxIter = 100, tol = 0.0001) {

	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("Family must be one of {'gaussian', 'binomial', 'poisson'} !")
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
	halpha 	= tbeta[1:p1]*scal_tx
	hbeta 	= tbeta[-c(1:p1)]*scal_x
	hdelta 	= (1+z%*%htheta>0)
	htheta 	= c(1, htheta)
	return(list(alpha = halpha, beta = hbeta, gamma = htheta, delta = hdelta))
}