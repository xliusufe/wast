exams <- function(family = "gaussian", method = "wast", tau = 0.5, B = 1000, K = 1000){
	if(family == "gaussian"){
		data(simulatedData_gaussian)
		pvals   = pvalglm(data = data_gaussian, family = family, method = method, B = B, K = K)
	}
	else if(family == "binomial"){
		data(simulatedData_binomial)
		pvals   = pvalglm(data = data_binomial, family = family, method = method, B = B, K = K)
	}
	else if(family == "poisson"){
		data(simulatedData_poisson)
		pvals   = pvalglm(data = data_poisson, family = family, method = method, B = B, K = K)
	}
	else if(family == "probit"){
		data(simulatedData_probit)
		pvals   = pval_probit(data = data_probit, method = method, B = B, K = K)
	}
	else if(family == "quantile"){
		data(simulatedData_quantile)
		pvals   = pval_quantile(data = data_quantile, method = method, tau = tau, B = B, K = K)
	}
	else if(family == "semiparam"){
		data(simulatedData_semiparam)
		pvals   = pval_semiparam(data = data_semiparam, method = method, B = B, K = K)
	}
	else{
		stop("Family must be one of {GLM with Gaussian family ('gaussian'), GLM with binomial family ('binomial'), GLM with Poisson family ('poisson'), Quantile regression ('quantile'), Probit regression ('probit'), and Semiparamtric models ('semiparam')} !")
	}

	return(pvals)
}