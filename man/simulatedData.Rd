\name{simulatedData}
\alias{simulatedData}
\docType{data}
\title{
Simulated data from generalized linear models
}
\description{
Simulated data from generalized linear models, including family 'gaussian' (\code{simulatedData_gaussian}), 
'binomial' (\code{simulatedData_binomial}), and 'poisson' (\code{simulatedData_poisson}). 
}
\usage{
data(simulatedData_gaussian)
}
\details{
We simulated data generated from generalized linear models
  \deqn{
  f(\boldsymbol{V}_i; \boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\gamma})=\exp\left\{ \frac{y_i\mu_i - c(\mu_i)}{a(\phi)}\right\}h(y_i),
  }
  where
  \deqn{
  \mu_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\gamma}\geq 0).
  }

\itemize{
\item Y: the response, an \eqn{n}-vector

\item X: the baseline variable with dimension \eqn{n\times p} 

\item Z: the grouping difference variable with dimension \eqn{n\times q} 

\item U: the grouping variable with dimension \eqn{n\times r} 
}

}

\references{
Liu, X. (2022). Subgroup detecting in generalized linear models. Manuscript.
}
\examples{
data(simulatedData_gaussian)

y <- data_gaussian$Y[1:5]
x <- dim(data_gaussian$X)
z <- dim(data_gaussian$Z)
u <- dim(data_gaussian$U)
}
\keyword{datasets}
