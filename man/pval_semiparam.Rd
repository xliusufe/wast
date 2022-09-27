\name{pval_semiparam}
\alias{pval_semiparam}
\title{
  P-value for subgroup test in semiparamtric models
}
\description{
  Provide p-value for subgroup test in semiparamtric models, including two methods 'wast' and 'sst'.
}

\usage{
  pval_semiparam(data, method = "wast", M = 1000, K = 2000, 
                    isBeta = FALSE, shape1 = 1, shape2 = 1)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{A} (treatment indicator), \eqn{X1} and \eqn{X2} (baseline variables), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}
  
  \item{method}{There are two methods, including the proposed 'wast' and 'sst'.}

  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst'.}

  \item{isBeta}{A bool value. The weight \eqn{w(\gamma)} is chosen to be Beta distribution if \code{isBeta=TRUE}, which can be used if the grouping difference variable is bounded in \eqn{[0,1]}. Default is \code{FALSE}.}

  \item{shape1}{The first parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{shape2}{The second parameter of Best distribution if \code{isBeta = TRUE}.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Semiparamtric models (see details in my paper Liu (2022))
  \deqn{
  Y_i = h(\boldsymbol{X}_{1i}) +A_i\boldsymbol{Z}_i^T\boldsymbol{\beta} \boldsymbol{1}(\boldsymbol{U}_i^T \boldsymbol{\theta}\geq0)+\epsilon_i,
  }
where \eqn{h(\boldsymbol{X}_1)} is an unknown baseline mean function for patients in treatment \eqn{A=0} which can be set a linear function \eqn{h(\boldsymbol{X}_1)  = X_1^T\boldsymbol{\alpha}_1}, 
and \eqn{P(A=1|\boldsymbol{X}_2)} can be modelled by a logistic regression model 
\deqn{
P(A=1|\boldsymbol{X}_2)=\pi(\boldsymbol{X}_2, \boldsymbol{\alpha}_2)=\frac{\exp\{\boldsymbol{X}_2^T\boldsymbol{\alpha}_2\}}{1+\exp\{\boldsymbol{X}_2^T\boldsymbol{\alpha}_2\}}.
}

The hypothesis test problem is 
\deqn{
  H_0: \boldsymbol{\beta}=\boldsymbol{0}\quad versus\quad H_1: \boldsymbol{\beta}\neq \boldsymbol{0}.
}
}



\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

LEE, S., SEO, M. H. and SHIN, Y. (2011). Testing for Threshold Effects in Regression Models. Journal of
the American Statistical Association 106, 220-231.

Liu, X. (2022). Change-plane testing in the generalized estimating equations. Manuscript.
}

\examples{ 
  data(simulatedData_semiparam)
  pvals <- pval_semiparam(data = data_semiparam, method = "wast")
  pvals

  pvals <- pval_semiparam(data = data_semiparam, method = "sst")
  pvals

}