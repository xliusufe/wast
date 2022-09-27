\name{pval_probit}
\alias{pval_probit}
\title{
  P-value for subgroup test in probit regression models
}
\description{
  Provide p-value for subgroup test in probit regression models, including two methods 'wast' and 'sst'.
}

\usage{
  pval_probit(data, method = "wast", M = 1000, K = 2000,
                  isBeta = FALSE, shape1 = 1, shape2 = 1)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

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
  Probit regression models
  \deqn{
  f(\boldsymbol{V}_i) = \Phi(h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta}))^{Y_i}+\Phi(-h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta}))^{1-Y_i},
  }
where \eqn{\Phi(\cdot)} is the cumulative distribution function of standard normal distribution,
and
\deqn{
h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta})
=\boldsymbol{X}_i\boldsymbol{\alpha}+\boldsymbol{Z}_i^T\boldsymbol{\beta}\boldsymbol{1}(\boldsymbol{U}_i^T\boldsymbol{\theta}\geq0).
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
  data(simulatedData_probit)
  pvals <- pval_probit(data = data_probit, method = "wast")
  pvals

  pvals <- pval_probit(data = data_probit, method = "sst")
  pvals

}