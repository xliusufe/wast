\name{pval_quantile}
\alias{pval_quantile}
\title{
  P-value for subgroup test in quantile regression models
}
\description{
  Provide p-value for subgroup test in quantile regression.
}

\usage{
  pval_quantile(data, method = "wast", tau = 0.5, M = 1000, K = 2000,
        isBeta = FALSE, shape1 = 1, shape2 = 1, N0 = 5000, MU = NULL, ZK = NULL)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{method}{There are two methods, including the proposed 'wast', 'wastapprox' and 'sst'.}

  \item{tau}{The given quantile \eqn{\tau}, a scale in the unit inteval. Default is \code{tau = 0.5}.}

  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values.}

  \item{isBeta}{A bool value. The weight \eqn{w(\gamma)} is chosen to be Beta distribution if \code{isBeta=TRUE}, which can be used if the grouping difference variable is bounded in \eqn{[0,1]}. Default is \code{FALSE}.}

  \item{shape1}{The first parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{shape2}{The second parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{N0}{An integer, the number of samples to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{N0 = 5000}. }

  \item{MU}{A vector with same length as \eqn{\boldsymbol{U}}, which is the mean of weight to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{MU = NULL}, in which \code{MU = runif(p3) - 0.5}. }

  \item{ZK}{A vector with length \code{N0}, which is normal sample to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{ZK = NULL}, in which \code{ZK	= rnorm(N0)}. }
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Quantile regression models
  \deqn{
  Q_{Y_i}(\tau|\boldsymbol{X}_i,\boldsymbol{Z}_i, \boldsymbol{U}_i)
  = \boldsymbol{X}_i^T\boldsymbol{\alpha}(\tau)+\boldsymbol{Z}_i^T\boldsymbol{\beta}(\tau)\boldsymbol{ 1}(\boldsymbol{U}_i^T\boldsymbol{\theta}(\tau)\geq0).
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
  data(simulatedData_quantile)
  pvals <- pval_quantile(data = data_quantile, tau = 0.5, method = "wast")
  pvals

  pvals <- pval_quantile(data = data_quantile, tau = 0.3, method = "wast")
  pvals

}