# SPSP: an R Package for Selecting the relevant predictors by Partitioning the Solution Paths of the Penalized Likelihood Approach

<!-- badges: start -->

[![](https://img.shields.io/cran/v/SPSP?logo=R)](https://cran.r-project.org/package=SPSP)
[![CRAN checks](https://badges.cranchecks.info/summary/SPSP.svg)](https://cran.r-project.org/web/checks/check_results_SPSP.html)
[![](https://cranlogs.r-pkg.org/badges/grand-total/SPSP?color=blue)](https://cranlogs.r-pkg.org/badges/grand-total/SPSP)
[![](https://cranlogs.r-pkg.org/badges/last-month/SPSP?color=green)](https://cranlogs.r-pkg.org/badges/last-month/SPSP?color=green)
[![](https://cranlogs.r-pkg.org/badges/last-week/SPSP?color=yellow)](https://cranlogs.r-pkg.org/badges/last-week/SPSP?color=yellow)

<!-- badges: end -->

Overview
--------

An implementation of the feature Selection procedure by Partitioning the entire Solution Paths
(namely SPSP) to identify the relevant features rather than using a single tuning parameter. 
By utilizing the entire solution paths, this procedure can obtain better selection accuracy than 
the commonly used approach of selecting only one tuning parameter based on existing criteria, 
cross-validation (CV), generalized CV, AIC, BIC, and EBIC (Liu, Y., & Wang, P. (2018) 
https://doi.org/10.1214/18-EJS1434). It is more stable and accurate (low false positive and false negative
rates) than other variable selection approaches. In addition, it can be flexibly coupled with 
the solution paths of Lasso, adaptive Lasso, SCAD, MCP, ridge regression, and other penalized estimators.

## Installation

The `SPSP` package is currently available on [SPSP CRAN](https://CRAN.R-project.org/package=SPSP).

### Install `SPSP` development version from GitHub (recommended)

``` r
# Install the development version from GitHub
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/SPSP")
```

### Install `SPSP` from the CRAN

``` r
# Install from CRAN
install.packages("SPSP")
```

## Example

The user-friendly function SPSP() conducts the selection by Partitioning the Solution
Paths (the SPSP procedure) to selects the relevant predictors. The user only needs 
to specify the independent variables matrix, response, family, and a penalized method
that can generate the solution paths, for example, Lasso, adaptive Lasso, SCAD, MCP, 
ridge regression. The embedded selection methods in this package can be called using
`fitfun.SP = lasso.glmnet`. Currently, six methods are included: `lasso.glmnet`,
`adalasso.glmnet`, `adalassoCV.glmnet`, `SCAD.ncvreg`, `MCP.ncvreg`, 
and `ridge.glmnet`.

The following example shows the R codes:

``` r
library(SPSP)
data(HihgDim)
library(glmnet)

x <- as.matrix(HighDim[,-1])
y <- HighDim[,1]

# SPSP + lasso
spsp_lasso_1 <- SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
                     init = 1, standardize = FALSE, intercept = FALSE)

head(spsp_lasso_1$nonzero)
head(spsp_lasso_1$beta_SPSP)

# SPSP + adalasso
spsp_adalasso_5 <- SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                        init = 5, standardize = T, intercept = FALSE)
                              
head(spsp_adalasso_5$nonzero)
head(spsp_adalasso_5$beta_SPSP)

# SPSP + SCAD
spsp_scad_5 <- SPSP(x = x, y = y, family = "gaussian", fitfun.SP = SCAD.ncvreg,
                    init = 5, standardize = T, intercept = FALSE)
                              
head(spsp_scad_5$nonzero)
head(spsp_scad_5$beta_SPSP)
```

References
----------

Liu, Y., & Wang, P. (2018). Selection by partitioning the solution paths. *Electronic Journal of Statistics*, 12(1), 1988-2017. <10.1214/18-EJS1434>

