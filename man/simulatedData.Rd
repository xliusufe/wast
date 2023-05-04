\name{simulatedData}
\alias{simulatedData}
\docType{data}
\title{
Simulated data from generalized linear models
}
\description{
Simulated data from the framework of general estimating equations, including model
\itemize{
\item 'Quantile regression' (\code{simulatedData_quantile}),
\item 'Probit regression' (\code{simulatedData_probit}),
\item 'Semiparamtric models' (\code{simulatedData_semiparam}),
\item Simulated data from generalized linear models, including family 'gaussian' (\code{simulatedData_gaussian}), 
'binomial' (\code{simulatedData_binomial}), and 'poisson' (\code{simulatedData_poisson}).
}
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
Liu, X. (2023). Subgroup testing in change-plane regression with high-dimensional grouping variables. Manuscript.
}
\examples{
data(simulatedData_gaussian)

y <- data_gaussian$Y[1:5]
x <- dim(data_gaussian$X)
z <- dim(data_gaussian$Z)
u <- dim(data_gaussian$U)

data(simulatedData_probit)
y <- data_probit$Y[1:5]
x <- dim(data_probit$X)
z <- dim(data_probit$Z)
u <- dim(data_probit$U)

data(simulatedData_semiparam)
y <- data_semiparam$Y[1:5]
x1 <- dim(data_semiparam$X1)
x2 <- dim(data_semiparam$X2)
z <- dim(data_semiparam$Z)
u <- dim(data_semiparam$U)
}
\keyword{datasets}