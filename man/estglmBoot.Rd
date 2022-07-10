\name{estglmBoot}
\alias{estglmBoot}
\title{
  Estimating standard deviation of parameters by bootstrap method in Generalized Linear Models with subgroups
}
\description{
  Provide estimators of standard deviation of coefficients by bootstrap method in generalized linear models with subgroups.
}

\usage{
  estglmBoot(data, family = "gaussian", h = NULL, smooth = "sigmoid", 
              weights = "exponential", B = 1000, maxIter = 100, tol = 0.0001)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{family}{Family for generalized linear models, including 'Gaussian', 'binomial', and 'Poisson'.}

  \item{h}{A numeric number, which is the bandwidth in the smooth function. Default is \code{h = log(n)/sqrt(n)}}

  \item{smooth}{The smooth function. Either "sigmoid" (the default), "pnorm", or "mixnorm", see details below.} 

  \item{weights}{The weights. Either "exponential" (the default), "norm", or "bernoulli", see details below.}

  \item{B}{An integer, the number of bootstrap samples. Default is \code{B = 1000}.}

  \item{maxIter}{An integer, the maximum number of iterations. Default is \code{maxIter = 100}.}

  \item{tol}{Convergence threshhold. Default is \code{tol = 0.0001}.}
}

\value{
  \item{alpha}{Estimator of the baseline parameter \eqn{\boldsymbol{\alpha}}.}
  \item{beta}{Estimator of the grouping difference parameter \eqn{\boldsymbol{\beta}}.}
  \item{gamma}{Estimator of the grouping parameter \eqn{\boldsymbol{\gamma}}.}
  \item{delta}{A vector with length \eqn{n}. Estimator of the indicator function \eqn{I(\boldsymbol{U}^T\boldsymbol{\gamma}\geq 0)}.}
  \item{std}{A vector with length \eqn{p+q+r-1}. The standard deviation (\code{sd}) of parameter \eqn{(\boldsymbol{\alpha}^T, \boldsymbol{\beta}^T,\boldsymbol{\gamma}_{-1}^T)^T}, 
  where \eqn{\boldsymbol{\gamma}_{-1}=(\gamma_2,\cdots,\gamma_r)^T}.}
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

  The smooth functioms: 
  \itemize{
  \item (a) sigmoid function ("sigmoid")
  \deqn{
  S(u)=1/(1+e^u);
  }
 
  \item (b) norm CDF ("pnorm")
   \deqn{
  S(u)=\Phi(u);
  } 

  \item (c) mixture of norm CDF and density ("mixnorm")
   \deqn{
  S(u)=\Phi(u) + u\phi(u),
  }
  }
  where \eqn{\Phi(u)} and \eqn{\phi(u)} are the CDF and density of starndard norm distribution, that is,
  \deqn{
  \Phi(u)=\int_{-\infty}^u\frac{1}{2\pi}\exp\left(-\frac{s^2}{2}\right)ds,
  } 
  and
  \deqn{
  \phi(u)=\frac{1}{2\pi}\exp\left(-\frac{u^2}{2}\right).
  }


  The weights from: 
    \describe{
    \item{(a)}{exponential distribution with unit rate parameter ("exponential");}
    \item{(b)}{normal distribution with unit mean and unit variance ("norm");}
    \item{(c)}{bernoulli distribution, of which value is 0 with probability 0.5 and 2 with probability 0.5.}
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
  fit <- estglmBoot(data = data_gaussian, family = "gaussian")
  fit$alpha
  fit$beta
  fit$theta
  fit$std
  data(simulatedData_binomial)
  fit <- estglmBoot(data = data_binomial, family = "binomial")
  fit$alpha
  fit$beta
  fit$theta
  fit$std

  data(simulatedData_poisson)
  fit <- estglmBoot(data = data_poisson, family = "poisson")
  fit$alpha
  fit$beta
  fit$theta
  fit$std
}