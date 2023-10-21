#' Selection by Partitioning the Solution Paths
#'
#' An implementation of the feature Selection procedure by Partitioning the entire Solution Paths
#' (namely SPSP) to identify the relevant features rather than using a single tuning parameter. 
#' By utilizing the entire solution paths, this procedure can obtain better selection accuracy than 
#' the commonly used approach of selecting only one tuning parameter based on existing criteria, 
#' cross-validation (CV), generalized CV, AIC, BIC, and EBIC (Liu, Y., & Wang, P. (2018)). It is 
#' more stable and accurate (low false positive and false negative rates) than other variable 
#' selection approaches. In addition, it can be flexibly coupled with the solution paths of Lasso, 
#' adaptive Lasso, ridge regression, and other penalized estimators.
#'
#' @details This package includes two main functions and several functions (fitfun.SP) to obtains
#' the solution paths. The \code{SPSP} function allows users to specify the penalized likelihood 
#' approaches that will generate the solution paths for the SPSP procedure. Then this function 
#' will automatically partitioning the entire solution paths. Its key idea is to classify variables
#' as relevant or irrelevant at each tuning parameter and then to select all of the variables 
#' which have been classified as relevant at least once. The \code{SPSP_step} purely apply the 
#' partitioning step that needs the solution paths as the input. In addition, there are several 
#' functions to obtain the solution paths. They can be used as an input of \code{fitfun.SP} argument.
#' 
#'
#' @docType package
#'
#' @author Xiaorui (Jeremy) Zhu, \email{zhuxiaorui1989@@gmail.com}, \cr Yang Liu, \email{yliu23@@fhcrc.org}, \cr
#' Peng Wang, \email{wangp9@@ucmail.uc.edu}
#' 
#' @references Liu, Y., & Wang, P. (2018). Selection by partitioning the solution paths. 
#' \emph{Electronic Journal of Statistics}, 12(1), 1988-2017. <10.1214/18-EJS1434>
#' 
#' @name SPSP-package
#' 
#' @useDynLib SPSP
#'
NULL

## Generic implementation of SPSP
SPSP <- function(x, ...) {
  UseMethod("SPSP", x)
}

#' Selection by partitioning the solution paths of Lasso, Adaptive Lasso, and Ridge penalized regression.
#' 
#' A user-friendly function to conduct the selection by Partitioning the Solution Paths (the SPSP algorithm). The
#' user only needs to specify the independent variables matrix, response, family, and \code{fitfun.SP}. 
#'
#' @param x A matrix with all independent variables, of dimension n by p; each row is an observation vector with p variables. 
#' @param y Response variable. Quantitative for \code{family="gaussian"} or \code{family="poisson"} (non-negative counts). 
#' For \code{family="binomial"} should be either a factor with two levels.
#' 
#' @param family Response type. Either a character string representing one of the built-in families,
#' or else a glm() family object.
#' @param fitfun.SP A function to obtain the solution paths for the SPSP algorithm. This function takes the arguments 
#' x, y, family as above, and additionally the standardize and intercept and others in \code{\link[glmnet]{glmnet}}, 
#' \code{\link[ncvreg]{ncvreg}}, or \code{\link[lars]{lars}}. The function fit the penalized models with lasso, adaptive lasso,
#' SCAD, or MCP penalty, or ridge regression to return the solution path of the corresponding penalized likelihood approach.
#' \describe{
#'   \item{\code{lasso.glmnet}}{lasso selection from \code{\link[glmnet]{glmnet}}.}
#'   \item{\code{adalasso.glmnet}}{adaptive lasso selection using the \code{lambda.1se} from cross-validation lasso method
#'   to obtain initial coefficients. It uses package \code{\link[glmnet]{glmnet}}.}
#'   \item{\code{adalassoCV.glmnet}}{adaptive lasso selection using the \code{lambda.1se} from cross-validation adaptive 
#'   lasso method to obtain initial coefficients. It uses package \code{\link[glmnet]{glmnet}}.}
#'   \item{\code{ridge.glmnet}}{use ridge regression to obtain the solution path.}
#'   \item{\code{SCAD.ncvreg}}{use SCAD-penalized regression model in \code{\link[ncvreg]{ncvreg}} to obtain the solution path.}
#'   \item{\code{MCP.ncvreg}}{use MCP-penalized regression model in \code{\link[ncvreg]{ncvreg}} to obtain the solution path.}
#'   \item{\code{lasso.lars}}{use lasso selection in \code{\link[lars]{lars}} to obtain the solution path.}
#' } 
#' @param standardize logical argument. Should conduct standardization before the estimation? Default is TRUE.
#' @param intercept logical. If x is a data.frame, this argument determines if the resulting model matrix should contain 
#' a separate intercept or not. Default is TRUE.
#' @param args.fitfun.SP A named list containing additional arguments that are passed to the fitting function;
#' see also argument \code{args.fitfun.SP} in do.call.
#' @param ... Additional optional arguments.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats glm.fit coef 
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom ncvreg ncvreg
#' 
#' @return An object of class \code{"SPSP"} is a list containing at least the following components:
#' \item{\code{beta_SPSP}}{the estimated coefficients of SPSP selected model;}
#' \item{\code{S0}}{the estimated relevant sets;}
#' \item{\code{nonzero}}{the selected covariates;}
#' \item{\code{zero}}{the covariates that are not selected;}
#' \item{\code{thres}}{the boundaries for abs(beta);}
#' \item{\code{R}}{the sorted adjacent distances;}
#' \item{\code{intercept}}{the estimated intercept when \code{intercept == T}.}
#' 
#' This object has attribute contains: 
#' 
#' \item{\code{mod.fit}}{the fitted penalized regression within the input function \code{fitfun.SP};}
#' \item{\code{family}}{the family of fitted object;}
#' \item{\code{fitfun.SP}}{the function to obtain the solution paths for the SPSP algorithm;}
#' \item{\code{args.fitfun.SP}}{a named list containing additional arguments for the function \code{fitfun.SP}.}
#' 
#' 
#' @export
#'
#' @examples
#' data(HighDim)
#' library(glmnet)
#' # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso and SPSP+AdaLasso:
#' data(HighDim)  
#' 
#' x <- as.matrix(HighDim[,-1])
#' y <- HighDim[,1]
#' 
#' spsp_lasso_1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
#'                            init = 1, standardize = FALSE, intercept = FALSE)
#' 
#' head(spsp_lasso_1$nonzero)
#' head(spsp_lasso_1$beta_SPSP)
#' 
#' spsp_adalasso_5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
#'                               init = 5, standardize = TRUE, intercept = FALSE)
#'                               
#' head(spsp_adalasso_5$nonzero)
#' head(spsp_adalasso_5$beta_SPSP)
#' 
#' 
#' 
SPSP <- 
  function(x, 
           y, 
           family = c("gaussian", "binomial"),
           fitfun.SP = adalasso.glmnet, 
           args.fitfun.SP = list(), 
           standardize = TRUE, 
           intercept = TRUE, 
           ...) {
  
  this.call <- match.call()
  family <- match.arg(family)
  
  glmnet_T <- inherits(fitfun.SP, "glmnet")
  
  fit_mod_SP <- do.call(fitfun.SP, c(list(x = x, y = y, family = family, 
                                          standardize = standardize, intercept = intercept),
                                     args.fitfun.SP))
  
  # Handle BETAs from different class of fitted models, beta may have different structure from different package.
  BETAs_tem <- getBETAs(fit_mod_SP)
  
  # Conduct the SPSP step that use the above solution path to obtain relevant predictors.
  SPSP_temp <- SPSP_step(x = x, y = y, BETA = BETAs_tem, standardize = standardize, intercept = intercept, ...)

  # Assign attributes and class
  attr(SPSP_temp, "mod.fit") <- fit_mod_SP
  attr(SPSP_temp, "family") <- family
  attr(SPSP_temp, "fitfun.SP") <- fitfun.SP
  attr(SPSP_temp, "args.fitfun.SP") <- args.fitfun.SP
  
  class(SPSP_temp) <- c("SPSP", class(fit_mod_SP), class(SPSP_temp))
  
  return(SPSP_temp)
}


#' The selection step with the input of the solution paths.
#' 
#' A function to select the relevant predictors by partitioning the solution paths (the SPSP algorithm) 
#' based on the user provided solution paths \code{BETA}. 
#'
#' @param x independent variables as a matrix, of dimension nobs x nvars; each row is an observation vector. 
#' @param y response variable. Quantitative for \code{family="gaussian"} or \code{family="poisson"} (non-negative counts). 
#' For \code{family="binomial"} should be either a factor with two levels.
#' 
#' @param family either a character string representing one of the built-in families, or else a glm() family object. 
#' @param BETA the solution paths obtained from a prespecified fitting step \code{fitfun.SP = lasso.glmnet} etc. It must be 
#' a p by k matrix, should be thicker and thicker, each column corresponds to a lambda, and lambda gets smaller and smaller.
#' It is the returned coefficient matrix \code{beta} from a \code{glmnet} object, or a \code{ncvreg} object.
#' @param standardize whether need standardization.
#' @param intercept logical. If x is a data.frame, this argument determines if the resulting model matrix should contain 
#' a separate intercept or not.
#' @param init initial coefficients, starting from init-th estimator of the solution paths. The default is 1. 
#' @param R sorted adjacent distances, default is NULL. Will be calculated inside.
#' @param ... Additional optional arguments.
#' 
#'
#' @examples 
#' data(HighDim)
#' library(glmnet)
#' 
#' x <- as.matrix(HighDim[,-1])
#' y <- HighDim[,1]
#' 
#' lasso_fit <- glmnet(x = x, y = y, alpha = 1, intercept = FALSE)
#' 
#' # SPSP+Lasso method
#' K <- dim(lasso_fit$beta)[2]
#' LBETA <- as.matrix(lasso_fit$beta)
#' 
#' spsp_lasso_1 <- SPSP_step(x = x, y = y, BETA = LBETA, 
#'                           init = 1, standardize = FALSE, intercept = FALSE)
#' head(spsp_lasso_1$nonzero)
#' head(spsp_lasso_1$beta_SPSP)
#' 
#' @return A list containing at least the following components:
#' \item{\code{beta_SPSP}}{the estimated coefficients of SPSP selected model;}
#' \item{\code{S0}}{the estimated relevant sets;}
#' \item{\code{nonzero}}{the selected covariates;}
#' \item{\code{zero}}{the covariates that are not selected;}
#' \item{\code{thres}}{the boundaries for abs(beta);}
#' \item{\code{R}}{the sorted adjacent distances;}
#' \item{\code{intercept}}{the estimated intercept when \code{intercept == T}.}
#' 
#' This object has attribute contains: 
#' 
#' \item{\code{mod.fit}}{the fitted penalized regression within the input function \code{fitfun.SP};}
#' \item{\code{family}}{the family of fitted object;}
#' \item{\code{fitfun.SP}}{the function to obtain the solution paths for the SPSP algorithm;}
#' \item{\code{args.fitfun.SP}}{a named list containing additional arguments for the function \code{fitfun.SP}.}
#' 
#' 
#' @export
#' 
SPSP_step <- 
  function(x, 
           y, 
           family = c("gaussian", "binomial"), 
           BETA,
           standardize=TRUE,
           intercept = TRUE, 
           init = 1, 
           R = NULL, 
           ...) {
  
  this.call <- match.call()
  family <- match.arg(family)
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (!is.logical(intercept)) 
    stop("Please specify intercept == TRUE/FALSE.")
  if (!is.logical(standardize)) 
    stop("Please specify stardardize == TRUE/FALSE.")
  
  
  if (standardize == TRUE) { # standardize == TRUE
    x <- scale(x, center = TRUE, scale = TRUE)
  }
  
  n <- dim(x)[1] # size
  p <- dim(x)[2] # dimension
  
  K <- dim(BETA)[2]
  BETA <- as.matrix(BETA[,K:1]) 
  # BETA for the following should be a p by k matrix;should be sparser and sparser
  # each column corresponds to a lambda, and lambda gets larger and larger
  
  K <- dim(BETA)[2] # the number of the tuning parameter lambda
  
  if (is.null(R)) {
    #### select the default the tuning parameter if it is not given
    gap0 <- sort(diff(c(0,sort(abs(BETA[,init])))), decreasing = TRUE)
    # sorted adjacent distances 
    R <- ifelse(gap0[2]!=0, as.numeric(gap0[1]/gap0[2]), as.numeric(gap0[1]))   
  }
  
  ### only consider the rows with at least one non-zero beta
  colsbeta <- colSums(BETA)
  K2 <- max(which(colsbeta!=0))
  
  thres <- c() ## boundaries for abs(beta)
  S0 <- list() # estimated relevant sets
  S0c <- list() # estimated irrelevant sets
  
  # the initial values; 
  thres_tmp <- max(abs(BETA[, 1]))
  S0_tmp <- which(abs(BETA[, 1]) > thres_tmp)
  S0c_tmp <- which(abs(BETA[, 1])  <=  thres_tmp)
  
  # loop for the update
  if (K2 >= 1) {
    for (k in 1:K2) {
      
      beta_abs <- abs(BETA[, k])
      beta_sort <- sort(beta_abs)
      
      # update the initials
      thres[k] <- ifelse(length(S0c_tmp) >= 1 ,max(beta_abs[S0c_tmp]),0)
      S0_tmp <- which(abs(BETA[, k]) > thres[k])
      S0c_tmp <- which(abs(BETA[, k]) <= thres[k])
      
      S0[[k]] <- S0_tmp
      S0c[[k]] <- S0c_tmp
      
      if (length(S0c[[k]]) >=1) {
        
        gap <- diff(c(0, beta_sort))
        
        # the distance between current relevant and irrelevant sets
        gap_10 <- ifelse(length(S0c_tmp) == p, 0, gap[length(S0c_tmp) + 1])
        
        # gap for the current irrelevant set: S0c_tmp
        gap_0 <- gap[1:length(S0c_tmp)]
        o1 <- which.max(gap_0)
        gap_01 <- max(gap_0)
        
        gap_02 <- ifelse(o1>1, max(gap[1:(o1-1)]), 0)
        
        if (gap_10  <=  R * gap_01 & gap_01 >= R * gap_02 ) {
          
          thres[k] <- ifelse(o1>1, beta_sort[o1-1], 0)
          
          S0_tmp <- which(abs(BETA[, k]) > thres[k])
          S0c_tmp <- which(abs(BETA[, k]) <= thres[k])
          
          S0[[k]] <- S0_tmp
          S0c[[k]] <- S0c_tmp
          
        }
      }
      
    }
  }
  
  index <- rep(1,p)
  
  for(i in 1:p) {
    if (all(abs(BETA[i, 1:K2]) <= thres)) index[i] <- 0
  }
  
  nz <- which(index == 1)
  z <- which(index == 0)
  beta_SPSP <- rep(0, p)
  
  if (intercept == TRUE) { # intercept == TRUE
    if (length(nz) >= 1) {
      xc <- x[, nz]
      xc1 <- cbind(1, xc)
      if (length(nz) < n) { # if p < n
        betac1 <- glm.fit(y = y, x = xc1, intercept = FALSE, family = family)$coefficients
      } else { # # if p >= n, use ridge reg to provide final coefficients
        # betac1 <- solve(t(xc1)%*%xc1+0.001*diag(length(nz)+1))%*%t(xc1)%*%y # Very slow
        # betac1 <- crossprod(solve(crossprod(xc1)+0.001*diag(length(nz)+1)), crossprod(xc1,y)) # Faster
        betac1 <- solve((crossprod(xc1) + 0.001 * diag(length(nz) + 1)), crossprod(xc1, y)) # Fastest
      }
      intercept_est <- betac1[1]
      betac <- betac1[-1]
      
      beta_SPSP[nz] <- betac  
    } else { # if no variable is selected and intercept == TRUE, run a model with intercept only
      intercept_est <- glm.fit(y = y, x = rep(1, n), intercept = FALSE, family = family)$coefficients
      betac <- NA
    } 
    
  } else { # intercept == FALSE
    intercept_est <- NA
    if (length(nz) >= 1) {
      xc <- x[, nz]
      
      if (length(nz) < n) {
        betac <- glm.fit(y = y, x = xc, intercept = FALSE, family = family)$coefficients
      } else { # change all as ridge reg to provide final coefficients
        # betac <- solve(t(xc)%*%xc+0.001*diag(length(nz)))%*%t(xc)%*%y # Very slow
        # betac <- crossprod(solve(crossprod(xc)+0.001*diag(length(nz))), crossprod(xc,y)) # Faster
        betac <- solve((crossprod(xc) + 0.001*diag(length(nz))), crossprod(xc, y)) # Fastest
      }
      
      beta_SPSP[nz] <- betac  
    }
    
    else if (!intercept) {
      stop("No variable is selected and no intercept.")
    }
  }  
  
  # Change the display and format of some results
  names(beta_SPSP) <- colnames(x)
  nz_name <- colnames(x)[nz]
  z_name <- colnames(x)[z]
  names(intercept_est) <- "(Intercept)"
  
  SPSP <- list(beta_SPSP = beta_SPSP, 
               S0 = S0,
               nonzero = nz_name, zero = z_name,
               thres = thres, R = R, 
               intercept = intercept_est)
  
  attr(SPSP, "thres") <- thres ## boundaries for abs(beta)
  attr(SPSP, "R") <- R ## sorted adjacent distances
  
  attr(SPSP, "y") <- y
  attr(SPSP, "x") <- x
  attr(SPSP, "init") <- init
  
  return(SPSP)
}


