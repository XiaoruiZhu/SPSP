################################################################################
# Generic function to transfer solution paths obtained from a prespecified fitting
# step \code{fitfun.SP = lasso.glmnet} to be a p by k matrix, should be thicker
# and thicker, each column corresponds to a lambda, and lambda gets smaller and smaller.
################################################################################

#' @keywords internal
getBETAs <- function(object, ...) {
  UseMethod("getBETAs")
}


#' @keywords internal
getBETAs.glmnet <- function(object, ...) {
  return(object$beta)
}



#' @importFrom Matrix Matrix
#' @keywords internal
getBETAs.ncvreg <- function(object, ...) {
  K <- dim(object$beta)[2]
  
  # Remove the intercept from beta matrix in ncvreg object
  BETAs <- Matrix(object$beta[-1,], sparse = TRUE)
  return(BETAs)
}
