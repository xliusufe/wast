exams <- function(family = "gaussian", method = "wast", M = 1000, K = 1000){
	if(family == "gaussian"){
		data(simulatedData_gaussian)
		pvals   = pval(data = data_gaussian, family = family, method = method, M=M, K = K)
	}
	else if(family == "binomial"){
		data(simulatedData_binomial)
		pvals   = pval(data = data_binomial, family = family, method = method, M=M, K = K)
	}
	else if(family == "poisson"){
		data(simulatedData_poisson)
		pvals   = pval(data = data_poisson, family = family, method = method, M=M, K = K)
	}
	else{
		stop("Family must be one of {'gaussian', 'binomial', 'poisson', 'gamma'} !")
	}

	return(pvals)
}