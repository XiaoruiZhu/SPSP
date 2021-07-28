test_that("Check SPSP_step() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  lasso_fit <- glmnet(x = x, y = y, alpha = 1, intercept = FALSE)
  
  # library(plotmo) # for plot_glmnet
  # plot_glmnet(lasso_fit, label=15, xvar ="lambda")
  # coef(lasso_fit, s = 0.5)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients
  
  
  # summary(test_lm_3)
  
  # SPSP+Lasso method
  K <- dim(lasso_fit$beta)[2]
  LBETA <- as.matrix(lasso_fit$beta)
  
  r_spsp1 <- SPSP::SPSP_step(x = x, y = y, BETA = LBETA, init = 1, standardize = FALSE, intercept = FALSE)
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP_step(x = x, y = y, BETA = LBETA, init = 5, standardize = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP_step(x = x, y = y, BETA = LBETA, init = 5, standardize = T, intercept = T)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  # r_spsp5_std$intercept
  
  # Comparison
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_3), est_beta_5_std[1:3])
  
  expect_equal(test_lm_int[1], r_spsp5_std$intercept)
})

test_that("Check SPSP() with lasso.glmnet() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients
  
  
  lasso_fit <- glmnet(x = x, y = y, alpha = 1, intercept = FALSE)
  # coef(lasso_fit)
  
  r_spsp1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
                        init = 1, standardize = FALSE, intercept = FALSE)
  
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  # head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
                        init = 5, standardize = FALSE, intercept = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                            init = 5, standardize = T, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  
  # Comparison
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_3), est_beta_5_std[1:3])
  
  expect_equal(test_lm_int[1], r_spsp5_std$intercept)
  
})

test_that("Check SPSP() with use data(Boston)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  library(MASS)
  # Boston Housing data from http://archive.ics.uci.edu/ml/datasets/Housing
  data(Boston, package="MASS")
  colnames(Boston) 
  summary(Boston)
  
  Boston2 <- subset(Boston, crim<=3.2)
  # Boston2 <- Boston
  # dim(Boston2)
  x <- as.matrix(Boston2[,-14]); y <- Boston2[,14] 
  summary(x)
  x[,"crim"] <- log(x[,"crim"])
  x[,"tax"] <- log(x[,"tax"])
  x[,"lstat"] <- log(x[,"lstat"])
  x[,"dis"] <- log(x[,"dis"])
  # x[,"rad"] <- log(x[,"rad"])
  y <- scale(y)
  x[,"age"] <- log(x[,"age"])
  
  for (i in 1:(ncol(x))){
    x[,i] <- scale(x[,i])
  }
  # summary(x)
  
  # round(cor(x), 2)
  # round(cor(Boston[,-14]), 2)
  # summary(y)
  
  expect_error(SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  expect_error(SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  expect_error(SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  expect_error(SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = ridge.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  err <- expect_error(SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = lasso.lars,
                          init = 5, standardize = T, intercept = FALSE))
  expect_equal(err$message, "The function lasso.lars() is currently under development.")
  
})

