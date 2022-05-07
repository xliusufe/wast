\name{pval}
\alias{pval}
\title{
  P-value for Subgroup Test in Generalized Linear Models
}
\description{
  Provide p-value for Subgroup Test in Generalized Linear Models, including three methods 'wast', 'sst', andd 'slrt'.
}

\usage{
  pval(data, family = "gaussian", method = 'wast', K = 2000, M=1000)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{family}{Family for generalized linear models, including 'Gaussian', 'binomial', and 'Poisson'.}
  
  \item{method}{There are there methods, including the proposed 'wast', 'sst', andd 'slrt'.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'slrt'.}
  
  \item{M}{An integer, the number of bootstrap samples.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Generalized linear models
  \eqn{
  $$
  \\
  f(\boldsymbol{V}_i; \boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\gamma})=\exp\left\{ \frac{y_i\mu_i - c(\mu_i)}{a(\phi)}\right\}h(y_i),\\
  $$
  $$
  \mu_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\gamma}\geq 0).
  \\
  $$
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
  pvals <- pval(data = data_gaussian, family = "gaussian")
  pvals

  data(simulatedData_binomial)
  pvals <- pval(data = data_binomial, family = "binomial")
  pvals

  data(simulatedData_poisson)
  pvals <- pval(data = data_poisson, family = "poisson")
  pvals
}