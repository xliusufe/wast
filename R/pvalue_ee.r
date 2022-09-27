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

EstTn_probit1 <- function(data, isBeta = 0, K = 2000L, M=1000L) {
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 50
	tol 	= 0.00001

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

	dims 	= c(n, p1, p2, p3, K, M, maxIter)
	fit 	<- .Call("_STT_Probit_bstr",
	 				as.integer(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(tol))
	pvals 	= fit$pvals
	alpha 	= fit$alpha
	theta 	= fit$theta

	return(pvals)
}

EstTn_probit0 <- function(data, isBeta = 0, shape1 = 1, shape2 = 1, K = 2000L, M=1000L) {
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 50
	tol 	= 0.0001
	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)
	params 	= c(shape1, shape2, tol)

	fit <- .Call("_Est_probit",
				as.integer(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(tol))
	alphahat = fit$coef
	mu 		= tx%*%alphahat
	resids 	= fit$residuals

	yb 		= matrix(0, n, M)
	for(k in 1:M){
		yb[,k] 	= rnorm(n)<mu
	}
	yb 	<- cbind(y, yb)
	fit <- .Call("_Probit_bstr",
				as.integer(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.integer(dims),
				as.numeric(params))
	teststat	= fit$Tn0
	teststat_p	= fit$Tn

	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

EstTn_semiparam1 <- function(data, K=2000L, M=1000L) {
	y 		= data$Y
	a 		= data$A
	n 		= length(y)
	tx1		= data$X1
	tx2		= data$X2
	x 		= data$Z
	z 		= data$U
	p11		= ifelse(is.null(ncol(tx1)) , 1, ncol(tx1))
	p12 	= ifelse(is.null(ncol(tx2)) , 1, ncol(tx2))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 50
	tol 	= 0.00001
	dims 	= c(n, p11, p12, p2, p3, K, M, maxIter)
	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}
	G   = matrix(rnorm(M*n), nrow = n, ncol = M)
	fit <- .Call("_SST_DoubleRobust_bstr",
				as.numeric(y),
				as.numeric(a),
				as.numeric(tx1),
				as.numeric(tx2),
				as.numeric(x),
				as.numeric(z),
				as.numeric(rtheta),
				as.numeric(G),
				as.integer(dims),
				as.numeric(tol))
	pvals 	= fit$pvals

	return(pvals)
}

EstTn_semiparam0 <- function(data, isBeta = 0, shape1 = 1, shape2 = 1, K=1000L, M=1000L) {
	# a~tx1 for the logistic model, that is E[a|tx1] = pi(tx1)
	# y~tx2 for linear model, that is, y = h(tx2) + eps
	y 		= data$Y
	a 		= data$A
	n 		= length(y)
	tx1		= data$X1
	tx2		= data$X2
	x 		= data$Z
	z 		= data$U

	p11		= ifelse(is.null(ncol(tx1)) , 1, ncol(tx1))
	p12 	= ifelse(is.null(ncol(tx2)) , 1, ncol(tx2))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 20
	tol 	= 0.0001
	params 	= c(shape1, shape2, tol)

	fit1	= .Call("_EST_LINEAR_H0",
					as.numeric(tx1),
					as.numeric(y),
					as.integer(n),
					as.integer(p11))
	alpha1 	= fit1$coef
	muhat1 	= tx1%*%alpha1
	resids1	= fit1$residuals

	dims 	= c(n, p11, p12, p2, p3, M, isBeta, maxIter)
	yb = matrix(0, n, M)
	for(k in 1:M){
		yb[,k] 	= muhat1 + resids1*rnorm(n)
	}
	yb = cbind(y, yb)
	fit <- .Call("_DoubleRobust_bstr",
				as.numeric(yb),
				as.numeric(a),
				as.numeric(tx1),
				as.numeric(tx1),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids1),
				as.integer(dims),
				as.numeric(params))

	teststat	= fit$Tn0
	teststat_p	= fit$Tn
	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

EstTn_quantile1 <- function(data, tau = 0.5, isBeta = 0, K = 20000L, M=1000L) {
	if(tau*(1-tau)<=0)	stop('tau is out of range!')
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))
	h 		= 2.0*n^(1/5)
	type 	= 2
	maxIter = 50
	tol 	= 0.00001

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}
	eb 		= matrix(sample(c(-2*tau, 2*(1-tau)), size=n*M, prob=c(tau,1-tau), replace=T), nrow = n, ncol = M)
	wb  	= matrix(rnorm(M*n), nrow = n, ncol = M)
	G 		= wb
	dims 	= c(n, p1, p2, p3, K, M, maxIter, type)
	params 	= c(tau, h, tol)
	fit 	<- .Call("_STT_Quantile_bstr",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params))
	pvals 	= fit$pvals
	alpha 	= fit$alpha
	theta 	= fit$theta

	return(pvals)
}

EstTn_quantile0 <- function(data, tau = 0.5, isBeta = 0, shape1 = 1, shape2 = 1, K = 2000L, M=1000L) {
	if(tau*(1-tau)<=0)	stop('tau is out of range!')
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 100
	tol 	= 0.00001
	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)
	params 	= c(tau, shape1, shape2, tol)

	fit <- .Call("_EST_QR",
				as.numeric(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(tau, tol)))
	alphahat = fit$coef
	muhat 	= tx%*%alphahat
	resids 	= fit$residuals

	yb = matrix(0, n, M)
	for(k in 1:M){
		wb 	= sample(c(-2*tau,2*(1-tau)), size=n, prob=c(tau,1-tau), replace=T)
		yb[,k] 	= muhat + abs(resids)*wb
	}
	yb 	<- cbind(y, yb)
	fit <- .Call("_Quantile_bstr",
				as.numeric(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.integer(dims),
				as.numeric(params))

	teststat	= fit$Tn0
	teststat_p	= fit$Tn
	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

EstTn_quantile_slrt <- function(data, tau = 0.5, K = 20000L, M=1000L, h = NULL) {
	if(tau*(1-tau)<=0)	stop('tau is out of range!')
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)) , 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	h 		= ifelse(is.null(ncol(h)) , 2.0*n^(1/5), h)
	maxIter = 50
	tol 	= 0.00001

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}
	wb  	= matrix(rnorm(M*n), nrow = n, ncol = M)
	G 		= wb
	dims 	= c(n, p1, p2, p3, K, M, maxIter)
	params 	= c(tau, tol, h)
	fit 	<- .Call("_Quantile_SLR",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params))
	pvals 	= fit$pvals[1]
	return(pvals)
}

pval_probit <- function(data, method = "wast", M = 1000, K = 2000, isBeta = FALSE, shape1 = 1, shape2 = 1){
	if(shape1<0 || shape2<0)
		stop("Two parameters of Beta distribution must not be negative!")
	isBeta = ifelse(isBeta, 1, 0)
	if(method=='sst'){
		pvals  	= EstTn_probit1(data, isBeta = isBeta, K = K, M=M)
	}
	else if(method=='wast'){
	   pvals  	= EstTn_probit0(data, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M=M)
	}
	else{
		stop("Method must be one of {'wast' and 'sst'} !")
	}
	return(pvals)
}

pval_semiparam <- function(data, method = "wast", M = 1000, K = 2000, isBeta = FALSE, shape1 = 1, shape2 = 1){
	if(shape1<0 || shape2<0)
		stop("Two parameters of Beta distribution must not be negative!")
	isBeta = ifelse(isBeta, 1, 0)
	if(method=='sst'){
		pvals  	= EstTn_semiparam1(data, K = K, M=M)
	}
	else if(method=='wast'){
	   pvals  	= EstTn_semiparam0(data, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M=M)
	}
	else{
		stop("Method must be one of {'wast' and 'sst'} !")
	}
	return(pvals)
}

pval_quantile <- function(data, method = "wast", tau = 0.5, M = 1000, K = 2000, isBeta = FALSE, shape1 = 1, shape2 = 1){
	if(shape1<0 || shape2<0)
		stop("Two parameters of Beta distribution must not be negative!")
	isBeta = ifelse(isBeta, 1, 0)
	if(method=='sst'){
		pvals  	= EstTn_quantile1(data, tau = tau, isBeta = isBeta, K = K, M=M)
	}
	else if(method=='wast'){
	   pvals  	= EstTn_quantile0(data, tau = tau, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M=M)
	}
	else if(method=='slrt'){
	   pvals  	= EstTn_quantile_slrt(data, tau = tau, K = K, M = M, h = NULL)
	}
	else{
		stop("Method must be one of {'wast', 'sst' and slrt} !")
	}

	return(pvals)
}
