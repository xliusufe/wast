\name{exams}
\alias{exams}
\title{
  Examples for Subgroup Test in Generalized Linear Models
}
\description{
  Examples for Family 'Gaussian', 'binomial', and 'Poisson'.
}

\usage{exams(family = "gaussian", method = "wast", tau = 0.5, M = 1000, K = 1000)
}

\arguments{
  \item{family}{Family for generalized linear models, including 'gaussian', 'binomial', 'poisson', 'probit', 'quantile', and 'semiparam'.}

  \item{method}{There are there methods, including the proposed 'wast', 'sst', and 'slrt'.}

  \item{tau}{The given quantile \eqn{\tau} when family is 'quantile', a scale in the unit inteval. Default is \code{tau = 0.5}.}

  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'slrt'.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
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
  pvals <- exams(family = "gaussian", method = "wast")
  pvals

  pvals <- exams(family = "binomial", method = "wast")
  pvals

  pvals <- exams(family = "poisson", method = "wast")
  pvals

  pvals <- exams(family = "quantile", method = "wast", tau = 0.5)
  pvals
}