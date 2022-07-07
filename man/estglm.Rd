\name{estglm}
\alias{estglm}
\title{
  Estimation in Generalized Linear Models with subgroups
}
\description{
  Provide estimators of coefficients in generalized linear models with subgroups.
}

\usage{
  estglm(data, family = "gaussian", h = NULL, maxIter = 100, tol = 0.0001)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{family}{Family for generalized linear models, including 'Gaussian', 'binomial', and 'Poisson'.}

  \item{h}{A numeric number, which is the bandwidth in the smooth function. Default is \code{h = log(n)/sqrt(n)}}

  \item{maxIter}{An integer, the maximum number of iterations. Default is \code{maxIter = 100}.}

  \item{tol}{Convergence threshhold. Default is \code{tol = 0.0001}.}
}

\value{
  \item{alpha}{Estimator of the baseline parameter \eqn{\boldsymbol{\alpha}}.}
  \item{beta}{Estimator of the grouping difference parameter \eqn{\boldsymbol{\beta}}.}
  \item{gamma}{Estimator of the grouping parameter \eqn{\boldsymbol{\gamma}}.}
  \item{delta}{A vector with length \eqn{n}. Estimator of the indicator function \eqn{I(\boldsymbol{U}^T\boldsymbol{\gamma}\geq 0)}.}
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
}

\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

Liu, X. (2022). Subgroup detecting in generalized linear models. Manuscript.
}

\examples{
  data(simulatedData_gaussian)
  fit <- estglm(data = data_gaussian, family = "gaussian")
  fit$alpha

  data(simulatedData_binomial)
  fit <- estglm(data = data_binomial, family = "binomial")
  fit$beta

  data(simulatedData_poisson)
  fit <- estglm(data = data_poisson, family = "poisson")
  fit$alpha
  fit$beta
}