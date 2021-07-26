test_that("SPSP Lasso works", {
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  test_lm <- lm(Y~., data = HighDim)
  test_lm_null <- lm(Y ~ 1, data = HighDim)
  # summary(test_lm)
  # model_step_f <- step(test_lm_null, scope=list(lower=test_lm_null, upper=test_lm), direction='forward')
  
  library(glmnet)
  lasso_fit <- glmnet(x = as.matrix(HighDim[,-1]), y = HighDim[,1], alpha = 1, intercept = FALSE)
  # summary(lasso_fit)
  library(plotmo) # for plot_glmnet
  plot_glmnet(lasso_fit, label=15, xvar ="lambda")
  
  coef(lasso_fit, s = 0.5)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  # summary(test_lm_3)
  
  # SPSP+Lasso method
  K <- dim(lasso_fit$beta)[2]
  LBETA <- as.matrix(lasso_fit$beta[,K:1])
  
  r_spsp1 <- SPSP::SPSP(x = HighDim[,-1], y = HighDim[,1], BETA = LBETA, init = 1, standardize = FALSE)
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  
  r_spsp5 <- SPSP::SPSP(x = HighDim[,-1], y = HighDim[,1], BETA = LBETA, init = 5, standardize = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = HighDim[,-1], y = HighDim[,1], BETA = LBETA, init = 5, standardize = T)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  # SPSP_names <- colnames(X)[which(est_beta_5 != 0)]
  
  
  # expect_equal(2 * 2, 4)
})
