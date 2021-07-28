#' A high dimensional dataset with n equals to 200 and p equals to 500. 
#' 
#' A dataset with 200 observations and 500 dimensions is generated from the following process:
#' linear regression model with only first three non-zero coefficients equal to 3, 2, and 1.5 respectively. 
#' The covariates are correlated with AR structure (rho=0.3). The error term is normally distributed with
#' zero mean and sd equals to 0.5.
#' 
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @name HighDim
#'
#' @usage
#' data(HighDim)
#'
#' @examples
#' 
#' # HighDim dataset is generated from the following process:
#' n <- 200; p <- 500; sigma <- 0.5
#' beta <- rep(0, p); nonzero <- c(1, 2, 3); zero <- setdiff(1:p, nonzero)
#' beta[nonzero] <- c(3, 2, 1.5)
#' Sigma <- 0.3^(abs(outer(1:p,1:p,"-")))
#' library(MASS)
#' X <- mvrnorm(n, rep(0,p), Sigma)
#' error <- rnorm(n, 0, sigma)
#' 
#' X <- apply(X, 2, scale) * sqrt(n)/sqrt(n-1)
#' error <- error - mean(error)
#' 
#' Y <- X %*% beta + error
#' HighDim <- data.frame(Y, X)
#' head(HighDim)
#' 
#' 
NULL
