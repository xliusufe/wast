\name{pvalglm}
\alias{pvalglm}
\title{
  P-value for Subgroup Test in Generalized Linear Models
}
\description{
  Provide p-value for subgroup test in generalized linear models, including three methods 'wast', 'sst', and 'slrt'.
}

\usage{
  pvalglm(data, family = "gaussian", method = 'wast', B = 1000, K = 2000,
          isBeta = FALSE, shape1 = 1, shape2 = 1, N0 = 5000, MU = NULL, ZK = NULL)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{family}{Family for generalized linear models, including 'gaussian', 'binomial', and 'poisson'.}

  \item{method}{There are there methods, including the proposed 'wast', 'wastapprox', 'sst', and 'slrt'.}

  \item{B}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'slrt'.}

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
  Generalized linear models
  \deqn{
  f(\boldsymbol{V}_i; \boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\gamma})=\exp\left\{ \frac{y_i\mu_i - c(\mu_i)}{a(\phi)}\right\}h(y_i),
  }
  where
  \deqn{
  \mu_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\gamma}\geq 0).
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

Liu, X. (2022). Change-plane testing in the generalized estimating equations. Manuscript.
}

\examples{ 
  data(simulatedData_gaussian)
  pvals <- pvalglm(data = data_gaussian, family = "gaussian")
  pvals

  data(simulatedData_binomial)
  pvals <- pvalglm(data = data_binomial, family = "binomial")
  pvals

  data(simulatedData_poisson)
  pvals <- pvalglm(data = data_poisson, family = "poisson")
  pvals
}