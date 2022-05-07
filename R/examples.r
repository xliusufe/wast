exams <- function(family = "gaussian", method = "wast", K = 1000, M = 1000){
	if(family == "gaussian"){
		data(simulatedData_gaussian)
		pvals   = pval(data = data_gaussian, family = family, method = method, K = K, M=M)
	}
	else if(family == "binomial"){
		data(simulatedData_binomial)
		pvals   = pval(data = data_binomial, family = family, method = method, K = K, M=M)
	}
	else if(family == "poisson"){
		data(simulatedData_poisson)
		pvals   = pval(data = data_poisson, family = family, method = method, K = K, M=M)
	}
	else{
		stop("Family must be one of {'gaussian', 'binomial', 'poisson', 'gamma'} !")
	}

	return(pvals)
}