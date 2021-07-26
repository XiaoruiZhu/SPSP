#' SPSP
#' 
#' Selecting the relevant predictors by Partitioning the Solution Paths (the SPSP algorithm)
#'
#' @param x independent variables as a matrix, of dimension nobs x nvars; each row is an observation vector. 
#' @param y response variable. Quantitative for \code{family="gaussian"} or \code{family="poisson"} (non-negative counts). 
#' For \code{family="binomial"} should be either a factor with two levels.
#' 
#' @param BETA a p by k matrix, should be sparser and sparser, each column corresponds to the estimated coefficients vector
#' for a given lambda, and lambda gets larger and larger.
#' @param family either a character string representing one of the built-in families, or else a glm() family object. 
#' @param init initial coefficients, starting from init-th estimator 
#' @param R sorted adjacent distances 
#' @param standardize  whether need standardization
#'
#' @importFrom Rcpp evalCpp
#' 
#' @return
#' @export
#'
#' @examples
#' data(HihgDim)
#' 
#' @useDynLib SPSP, .registration=TRUE
SPSP <- function(x, y, BETA, family = gaussian, init=1, R=NULL, standardize=FALSE){
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  n <- dim(x)[1] # size
  p <- dim(x)[2] # dimension
  BETA <- as.matrix(BETA) 
  # BETA should be a p by k matrix;should be sparser and sparser
  # each column corrspends to a lambda, and lambda gets larger and larger
  K <- dim(BETA)[2] # the number of the tuning parameter lambda
  
  if(is.null(R)){
    #### select the default the tuning parameter if it is not given
    gap0 <- sort(diff(c(0,sort(abs(BETA[,init])))),decreasing = TRUE)
    # sorted adjacent distances 
    R <- ifelse(gap0[2]!=0,as.numeric(gap0[1]/gap0[2]),as.numeric(gap0[1]))   
  }
  
  ### only consider the rows with at least one non-zero beta
  colsbeta <- colSums(BETA)
  K2 <- max(which(colsbeta!=0))
  
  thres <- c() ## boundaries for abs(beta)
  S0 <- list() # estimated relevant sets
  S0c <- list() # estimated irrelevant sets
  
  # the initial values; 
  thres_tmp <- max(abs(BETA[,1]))
  S0_tmp <- which(abs(BETA[,1])>thres_tmp)
  S0c_tmp <- which(abs(BETA[,1])<=thres_tmp)
  
  # loop for the update
  if(K2>=1){
    for (k in 1:K2){
      
      beta_abs <- abs(BETA[,k])
      beta_sort <- sort(beta_abs)
      
      # update the initials
      thres[k] <- ifelse(length(S0c_tmp)>=1,max(beta_abs[S0c_tmp]),0)
      S0_tmp <- which(abs(BETA[,k])>thres[k])
      S0c_tmp <- which(abs(BETA[,k])<=thres[k])
      
      S0[[k]] <- S0_tmp
      S0c[[k]] <- S0c_tmp
      
      if(length(S0c[[k]])>=1){
        
        gap <- diff(c(0,beta_sort))
        
        # the distance between current relevant and irrelevant sets
        gap_10 <- ifelse(length(S0c_tmp)==p,0,gap[length(S0c_tmp)+1])
        
        # gap for the current irrelevant set: S0c_tmp
        gap_0 <- gap[1:length(S0c_tmp)]
        o1 <- which.max(gap_0)
        gap_01 <- max(gap_0)
        
        gap_02 <- ifelse(o1>1,max(gap[1:(o1-1)]),0)
        
        if(gap_10 <= R * gap_01 & gap_01 >= R*gap_02 ){
          
          thres[k] <- ifelse(o1>1,beta_sort[o1-1],0)
          
          S0_tmp <- which(abs(BETA[,k])>thres[k])
          S0c_tmp <- which(abs(BETA[,k])<=thres[k])
          
          S0[[k]] <- S0_tmp
          S0c[[k]] <- S0c_tmp
          
        }
      }
      
      
    }
  }
  
  
  index <- rep(1,p)
  
  for(i in 1:p){
    if(all(abs(BETA[i,1:K2])<=thres)) index[i] <- 0
  }
  
  
  nz <- which(index==1)
  z <- which(index==0)
  beta_SPSP <- rep(0,p)
  
  if(standardize==FALSE){
    intercept <- 0
    if(length(nz)>=1){
      xc <- x[,nz]
      
      if(length(nz)<=n){
        # betac <- .lm.fit(y~xc-1)$coefficients
        betac <- glm.fit(y = y, x = xc, intercept = FALSE, family = family)$coefficients
      }else{ # change all as ridge reg to provide final coefficients
        # betac <- solve(t(xc)%*%xc+0.001*diag(length(nz)))%*%t(xc)%*%y # Very slow
        betac <- solve((crossprod(xc) + 0.001*diag(length(nz))), crossprod(xc, y))
        # betac <- crossprod(solve(crossprod(xc)+0.001*diag(length(nz))), crossprod(xc,y)) # Faster
      }
      
      beta_SPSP[nz] <- betac  
    }
  }else if(standardize==TRUE){
    if(length(nz)>=1){
      xc <- x[,nz]
      xc1 <- cbind(1,xc)
      if(length(nz)<n){
        # betac1 <- (.lm.fit(y~xc1-1)$coefficients)
        betac1 <- glm.fit(y = y, x = xc1, intercept = FALSE, family = family)$coefficients
        
        intercept <- betac1[1]
        betac <- betac1[-1]
      }else{ # change all as ridge reg to provide final coefficients
        # betac1 <- solve(t(xc1)%*%xc1+0.001*diag(length(nz)+1))%*%t(xc1)%*%y # Very slow
        # betac1 <- crossprod(solve(crossprod(xc1)+0.001*diag(length(nz)+1)), crossprod(xc1,y)) # Faster
        betac1 <- solve((crossprod(xc1) + 0.001*diag(length(nz))), crossprod(xc1, y)) # Fastest
        intercept <- betac1[1]
        betac <- betac1[-1]
      }
      
      beta_SPSP[nz] <- betac  
    }else{
      # intercept <- .lm.fit(y~1)$coefficients[1]
      intercept <- glm.fit(y = y, x = rep(1,n), intercept = TRUE, family = family)$coefficients
    }
  }else{
    print("Please specify stardardize==TRUE/FALSE.")
  }
  
  
  list(beta_SPSP = beta_SPSP, S0 = S0, thres = thres,
       nonzero=nz, zero=z, R=R, intercept=as.numeric(intercept))
}


