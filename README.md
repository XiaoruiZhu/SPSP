# SPSP: an R Package for Selecting the relevant predictors by Partitioning the Solution Paths of the Penalized Likelihood Approach

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/SPSP)](https://www.r-pkg.org/badges/version/SPSP)
[![CRAN checks](https://cranchecks.info/badges/summary/SPSP)](https://cran.r-project.org/web/checks/check_results_SPSP.html)
[![](https://cranlogs.r-pkg.org/badges/grand-total/SPSP?color=blue)](https://cranlogs.r-pkg.org/badges/grand-total/SPSP)
[![](https://cranlogs.r-pkg.org/badges/last-month/SPSP?color=green)](https://cranlogs.r-pkg.org/badges/last-month/SPSP?color=green)
[![](https://cranlogs.r-pkg.org/badges/last-week/SPSP?color=yellow)](https://cranlogs.r-pkg.org/badges/last-week/SPSP?color=yellow)
[![](https://api.travis-ci.com/XiaoruiZhu/SPSP.svg?branch=master)](https://api.travis-ci.com/XiaoruiZhu/SPSP.svg)

<!-- badges: end -->

Overview
--------

An implementation of the feature Selection procedure by Partitioning the entire Solution Paths
(namely SPSP) to identify the relevant features rather than using a single tuning parameter. 
By utilizing the entire solution paths, this procedure can obtain better selection accuracy than 
the commonly used approach of selecting only one tuning parameter based on existing criteria, 
cross-validation (CV), generalized CV, AIC, BIC, and EBIC (Liu, Y., & Wang, P. (2018)). It is 
more stable and accurate (low false positive and false negative rates) than other variable 
selection approaches. In addition, it can be flexibly coupled with the solution paths of Lasso, 
adaptive Lasso, ridge regression,  and other penalized estimators.

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

The following example shows the R code for 

``` r
library(SPSP)
data(HihgDim)
library(glmnet)

# Use the user-friendly function SPSP() to conduct the selection by Partitioning the 
# Solution Paths (the SPSP algorithm). The user only needs to specify the independent 
# variables matrix, response, family, and \code{fitfun.SP = lasso.glmnet}. 

x <- as.matrix(HighDim[,-1])
y <- HighDim[,1]

spsp_lasso_1 <- SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
                     init = 1, standardize = FALSE, intercept = FALSE)

head(spsp_lasso_1$nonzero)
head(spsp_lasso_1$beta_SPSP)

spsp_adalasso_5 <- SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                        init = 5, standardize = T, intercept = FALSE)
                              
head(spsp_adalasso_5$nonzero)
head(spsp_adalasso_5$beta_SPSP)

# Use the function SPSP_step() to select the relevant predictors by partitioning the 
# solution paths based on the user provided solution paths \code{BETA}. 

lasso_fit <- glmnet(x = x, y = y, alpha = 1, intercept = FALSE)

# SPSP+Lasso method
K <- dim(lasso_fit$beta)[2]
LBETA <- as.matrix(lasso_fit$beta)
spsp_lasso_1 <- SPSP_step(x = x, y = y, BETA = LBETA, 
                          init = 1, standardize = FALSE, intercept = FALSE)

head(spsp_lasso_1$nonzero)
head(spsp_lasso_1$beta_SPSP)
```

References
----------

Liu, Y., & Wang, P. (2018). Selection by partitioning the solution paths. *Electronic Journal of Statistics*, 12(1), 1988-2017. <10.1214/18-EJS1434>

