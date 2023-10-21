################################################################################
## fitfun.SP functions
##
## a function to obtain the solution paths for the SPSP algorithm. These 
## functions need to take arguments x, y, family, standardize, and intercept 
## and return a glmnet object that has the solution paths stored. The solution 
## paths will be used as ingredients for the SPSP algorithm.
##
################################################################################



#' @title Four Fitting-Functions that can be used as an input of \code{fitfun.SP} argument 
#' to obtain the solution paths for the SPSP algorithm. The users can also customize a
#' function to generate the solution paths. As long as the customized function take 
#' arguments x, y, family, standardize, and intercept, and return an object of class 
#' \code{glmnet}, \code{lars} (or \code{SCAD}, \code{MCP} in the future).
#'
#' @name Fitting-Functions
NULL



#' \code{lasso.glmnet} uses lasso selection from \code{\link[glmnet]{glmnet}}.
#'
#' @param x a matrix of the independent variables. The dimensions are (nobs) and (nvars); each row is an observation vector. 
#' @param y Response variable. Quantitative for \code{family="gaussian"} or \code{family="poisson"} (non-negative counts). 
#' For \code{family="binomial"} should be either a factor with two levels.
#' 
#' @param family Response type. Either a character string representing one of the built-in families,
#' or else a glm() family object.
#' @param standardize logical argument. Should conduct standardization before the estimation? Default is TRUE.
#' @param intercept logical. If x is a data.frame, this argument determines if the resulting model matrix should contain 
#' a separate intercept or not. Default is TRUE.
#' @param ... Additional optional arguments.
#' 
#' @return An object of class \code{"glmnet"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @rdname Fitting-Functions
#' 
#' @importFrom glmnet glmnet
#' 
#' @export
#'
lasso.glmnet <- function(x, 
                         y, 
                         family,
                         standardize, 
                         intercept, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  fit_sp <- glmnet(x = x, y = y, family = family, 
                   alpha=1, standardize = standardize, intercept=intercept, ...) 
  return(fit_sp)
}


#' \code{adalasso.glmnet} the function to conduct the adaptive lasso selection using the \code{lambda.1se} 
#' from cross-validation lasso method to obtain initial coefficients. It uses package \code{\link[glmnet]{glmnet}}.
#' @rdname Fitting-Functions
#' 
#' @return An object of class \code{"glmnet"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @importFrom glmnet glmnet cv.glmnet
#' 
#' @export
#'
adalasso.glmnet <- function(x, 
                            y, 
                            family,
                            standardize, 
                            intercept, ...) {
  # Extract all optional arguments for further use. 
  args <- list(...)
  for (name in names(args) ) {
    assign(name, args[[name]])
  }
  
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " is needed but not available")
  
  n <- dim(x)[1]; p <- dim(x)[2]
  
  if (n<p) {
    # use CV for the ridge regression and "temp$lambda.1se" to obtain initial coefficients
    cv_fl1 <- cv.glmnet(x = x, y = y, family = family,
                        alpha = 0,
                        intercept = intercept, standardize = standardize, ...)
    
    first_coef <- coef(cv_fl1, s = cv_fl1$lambda.1se)
    first_beta <- first_coef[-1]
    
  } else {
    if (intercept) {
      x1 <- as.matrix(cbind(1, x))
      # first_coef <- solve(t(x)%*%x)%*%t(x)%*%y
      first_coef <- crossprod(solve(crossprod(x1)), crossprod(x1,y))
      first_beta <- first_coef[-1]
    } else {
      x1 <- as.matrix(x)
      # first_coef <- solve(t(x)%*%x)%*%t(x)%*%y
      first_coef <- crossprod(solve(crossprod(x1)), crossprod(x1,y))
      first_beta <- first_coef
    }
  }
  
  # If the penalty.factor are given, then adjust the penalties.
  
  if (exists("penalty.factor")) {
    penalty.factor <- abs(first_beta + 1/sqrt(n))^(-1) * penalty.factor
  } else {
    penalty.factor <- abs(first_beta + 1/sqrt(n))^(-1)
  }
  
  fit_sp <- glmnet(x = x, y = y, #penalty.factor = penalty.factor, 
                   family = family, alpha=1, 
                   intercept = intercept, standardize = standardize, ...)
  
  return(fit_sp)
}

#' \code{adalassoCV.glmnet} adaptive lasso selection using the \code{lambda.1se} from cross-validation adaptive 
#'   lasso method to obtain initial coefficients. It uses package \code{\link[glmnet]{glmnet}}.
#' @rdname Fitting-Functions
#' 
#' @return An object of class \code{"glmnet"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats coef
#' 
#' @export
#'
adalassoCV.glmnet <- function(x, 
                              y, 
                              family,
                              standardize, 
                              intercept, ...) {
  # Extract all optional arguments for further use. 
  args <- list(...)
  for (name in names(args) ) {
    assign(name, args[[name]])
  }
  
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  n <- dim(x)[1]; p <- dim(x)[2]
  
  if (n<p) {
    fl1 <- glmnet(x = x, y = y, family = family,
                  alpha=1, 
                  intercept = intercept,
                  standardize = standardize, ...)
    # use lasso from cv.glmnet to find initial lambda
    cv_fl1 <- cv.glmnet(x = x, y = y, family = family,
                        alpha=1, intercept = intercept, standardize = standardize, ...)
    lambda_tem <- cv_fl1$lambda.1se # use this as prespecified lambda
    first_coef <- fl1$beta[,which.min(abs(fl1$lambda-lambda_tem))]
    
  } else {
    # first_coef <- solve(t(x)%*%x)%*%t(x)%*%y
    first_coef <- crossprod(solve(crossprod(x)), crossprod(x,y))
  }
  
  # If the penalty.factor are given, then adjust the penalties.
  if (exists("penalty.factor")) {
    penalty.factor <- abs(first_coef + 1/sqrt(n))^(-1) * penalty.factor
  } else {
    penalty.factor <- abs(first_coef + 1/sqrt(n))^(-1)
  }
  
  fit_cv_ada <- cv.glmnet(x = x, y = y, # penalty.factor = penalty.factor, 
                   family = family, alpha=1, intercept = intercept, standardize = standardize, ...) 
  
  first_ada_coef <- coef(fit_cv_ada, s=fit_cv_ada$lambda.1se)
  first_ada_beta <- first_ada_coef[-1]
  
  
  # If the penalty.factor are given, then adjust the penalties.
  if (exists("penalty.factor")) {
    penalty.factor <- abs(first_ada_beta + 1/sqrt(n))^(-1) * penalty.factor
  } else {
    penalty.factor <- abs(first_ada_beta + 1/sqrt(n))^(-1)
  }
  
  fit_sp <- glmnet(x = x, y = y, # penalty.factor = penalty.factor, 
                   family = family, alpha=1, 
                   intercept = intercept, standardize = standardize, ...)
  
  # Store the first adalasso coefficient obtained from the adalasso.CV procedure. 
  fit_sp$first_ada_coef <- first_ada_coef
  
  return(fit_sp)
}

#' \code{ridge.glmnet} uses ridge regression to obtain the solution path.
#'
#' @return An object of class \code{"glmnet"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @rdname Fitting-Functions
#' @export
#'
ridge.glmnet <- function(x, 
                         y, 
                         family,
                         standardize, 
                         intercept, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " is needed but not available")
  
  fit_sp <- glmnet(x = x, y = y, family = family, 
                   alpha = 0, standardize = standardize, intercept=intercept, ...) 
  return(fit_sp)
}

#' \code{lasso.lars} uses lasso selection in \code{\link[lars]{lars}} to obtain the solution path.
#' 
#' @rdname Fitting-Functions
#' 
#' @return An object of class \code{"lars"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @importFrom lars lars
#' 
#' @export
#'
lasso.lars <- function(x, 
                       y, 
                       family,
                       standardize, 
                       intercept, ...) {
  if (!requireNamespace("lars", quietly = TRUE)) 
    stop("Package ", sQuote("lars"), " needed but not available")
  
  stop("The function lasso.lars() is currently under development.")
  
  # return(fit_sp)
}

#' \code{SCAD.ncvreg} uses SCAD penalty from \code{\link[ncvreg]{ncvreg}} for fitting regularization paths.
#' 
#' @rdname Fitting-Functions
#' 
#' @return An object of class \code{"ncvreg"} is returned to provide SCAD penalty solution paths for the SPSP algorithm. 
#' 
#' @importFrom ncvreg ncvreg
#' 
#' @export
#'
SCAD.ncvreg <- function(x, 
                       y, 
                       family,
                       standardize, 
                       intercept, ...) {
  if (!requireNamespace("ncvreg", quietly = TRUE)) 
    stop("Package ", sQuote("ncvreg"), " is needed but not available")
  
  # If the penalty.factor are given, use the penalties.
  
  fit_sp <- ncvreg(X = x,
                   y = y,
                   family = family,
                   alpha = 1,
                   penalty = "SCAD",
                   ...)
  
  return(fit_sp)
}

#' \code{MCP.ncvreg} uses MCP penalty from \code{\link[ncvreg]{ncvreg}} for fitting regularization paths.
#' 
#' @rdname Fitting-Functions
#' 
#' @return An object of class \code{"ncvreg"} is returned to provide solution paths for the SPSP algorithm. 
#' 
#' @importFrom ncvreg ncvreg
#' 
#' @export
#'
MCP.ncvreg <- function(x, 
                        y, 
                        family,
                        standardize, 
                        intercept, ...) {
  if (!requireNamespace("ncvreg", quietly = TRUE)) 
    stop("Package ", sQuote("ncvreg"), " is needed but not available")
  
  # If the penalty.factor are given, use the penalties.
  
  fit_sp <- ncvreg(X = x,
                   y = y,
                   family = family,
                   alpha = 1,
                   penalty = "MCP",
                   ...)
  
  return(fit_sp)
}