# wast
R package "wast" for calculating p-value of the test statistic for subgroup detecting in generalized linear models. In the paper Liu (2022), we consider hypothesis test of coefficients in the generalized linear models (GLM) to detect the existence of the subgroups, which can serve as the optimal individualized treatment recommendation in practice. Test we consider in this paper is one of the class of test problems when a part of parameters is not identifiable under the null. We propose a novel U-like statistic by taking the weighted average over the grouping parameter's space. The proposed test statistic not only improves significantly the power but also is computationally efficient.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/wast")

# Usage

   - [x] [wast-manual.pdf](https://github.com/xliusufe/wast/blob/master/inst/wast-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(wast)

    data(simulatedData_gaussian)
    pvals <- pval(data = data_gaussian, family = "gaussian")
    pvals

    data(simulatedData_binomial)
    pvals <- pval(data = data_binomial, family = "binomial")
    pvals

    data(simulatedData_poisson)
    pvals <- pval(data = data_poisson, family = "poisson")
    pvals

# References
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Davies, R. B. (1977). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 64(2):247-254.

Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 74(1):33-43.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 70(2):291-311.

Liu, X. (2022). Subgroup detecting in generalized linear models. Manuscript.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
