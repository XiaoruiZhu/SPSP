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
  # expect_equal(coef(test_lm_3), est_beta_5_std[1:3])
  
  expect_equal(test_lm_int[1], r_spsp5_std$intercept)
})

test_that("Check SPSP() with lasso.glmnet() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  xstd <- scale(HighDim[,-1], center = TRUE, scale = TRUE)
  
  HighDim2 <- data.frame(Y = y, xstd)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  test_lm_std_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim2)
  
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
                            init = 5, standardize = TRUE, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  r_spsp5_std_int <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                            init = 5, standardize = TRUE, intercept = TRUE)
  est_beta_5_std_int <- r_spsp5_std_int$beta_SPSP
  # head(est_beta_5_std_int)
  # r_spsp5_std_int$intercept
  
  # Check if coefficients are same as ols
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std_int[1:3])
  
  # Check if intercept is correct
  expect_true(is.na(r_spsp1$intercept))
  expect_true(is.na(r_spsp5$intercept))
  expect_true(is.na(r_spsp5_std$intercept))
  expect_true(!is.na(r_spsp5_std_int$intercept))
})

test_that("Check SPSP() with adalasso.glmnet() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  xstd <- scale(HighDim[,-1], center = TRUE, scale = TRUE)
  
  HighDim2 <- data.frame(Y = y, xstd)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  test_lm_std_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim2)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients
  
  r_spsp1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                        init = 1, standardize = FALSE, intercept = FALSE)
  
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  # head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                        init = 5, standardize = FALSE, intercept = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                            init = 5, standardize = TRUE, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  r_spsp5_std_int <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalasso.glmnet,
                                init = 5, standardize = TRUE, intercept = TRUE)
  est_beta_5_std_int <- r_spsp5_std_int$beta_SPSP
  # head(est_beta_5_std_int)
  # r_spsp5_std_int$intercept
  
  # Check if coefficients are same as ols
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std_int[1:3])
  
  # Check if intercept is correct
  expect_true(is.na(r_spsp1$intercept))
  expect_true(is.na(r_spsp5$intercept))
  expect_true(is.na(r_spsp5_std$intercept))
  expect_true(!is.na(r_spsp5_std_int$intercept))
})

test_that("Check SPSP() with adalassoCV.glmnet() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  xstd <- scale(HighDim[,-1], center = TRUE, scale = TRUE)
  
  HighDim2 <- data.frame(Y = y, xstd)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  test_lm_std_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim2)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients
  
  r_spsp1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                        init = 1, standardize = FALSE, intercept = FALSE)
  
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  # head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                        init = 5, standardize = FALSE, intercept = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                            init = 5, standardize = TRUE, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  r_spsp5_std_int <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                                init = 5, standardize = TRUE, intercept = TRUE)
  est_beta_5_std_int <- r_spsp5_std_int$beta_SPSP
  # head(est_beta_5_std_int)
  # r_spsp5_std_int$intercept
  
  # Check if coefficients are same as ols
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std_int[1:3])
  
  # Check if intercept is correct
  expect_true(is.na(r_spsp1$intercept))
  expect_true(is.na(r_spsp5$intercept))
  expect_true(is.na(r_spsp5_std$intercept))
  expect_true(!is.na(r_spsp5_std_int$intercept))
})

test_that("Check SPSP() with SCAD.ncvreg() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(ncvreg)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  xstd <- scale(HighDim[,-1], center = TRUE, scale = TRUE)
  
  HighDim2 <- data.frame(Y = y, xstd)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  test_lm_std_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim2)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients

  test <- SCAD.ncvreg(x = x, 
                          y = y, 
                          family = "gaussian",
                          standardize = FALSE, 
                          intercept = FALSE)
  # test$beta
  
  r_spsp1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = SCAD.ncvreg,
                        init = 1, standardize = FALSE, intercept = FALSE)
  
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  # head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = SCAD.ncvreg,
                        init = 5, standardize = FALSE, intercept = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = SCAD.ncvreg,
                            init = 5, standardize = TRUE, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  r_spsp5_std_int <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = SCAD.ncvreg,
                                init = 5, standardize = TRUE, intercept = TRUE)
  est_beta_5_std_int <- r_spsp5_std_int$beta_SPSP
  # head(est_beta_5_std_int)
  # r_spsp5_std_int$intercept
  
  # Check if coefficients are same as ols
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std_int[1:3])
  
  # Check if intercept is correct
  expect_true(is.na(r_spsp1$intercept))
  expect_true(is.na(r_spsp5$intercept))
  expect_true(is.na(r_spsp5_std$intercept))
  expect_true(!is.na(r_spsp5_std_int$intercept))
})

test_that("Check SPSP() with MCP.ncvreg() use data(HighDim)", {
  skip_if_not_installed("glmnet")
  
  library(ncvreg)
  
  # Use the high dimensional dataset (data(HighDim)) to test SPSP+Lasso:
  data(HighDim)  
  
  x <- as.matrix(HighDim[,-1])
  y <- HighDim[,1]
  
  xstd <- scale(HighDim[,-1], center = TRUE, scale = TRUE)
  
  HighDim2 <- data.frame(Y = y, xstd)
  
  test_lm_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim)
  
  test_lm_std_3 <- lm(Y ~ X1 + X2 + X3 + 0, data = HighDim2)
  
  # When nonzero == 0, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = rep(1, length(y)), intercept = TRUE, family = gaussian())$coefficients
  
  # When nonzero >= 1, run the following to obtain intercept
  # test_lm_int <- glm.fit(y = y, x = cbind(1, x[,1:3]), intercept = FALSE, family = gaussian())$coefficients
  test_lm_int <- glm(y ~ X1 + X2 + X3 + 1, data = HighDim, family = gaussian())$coefficients
  
  test <- MCP.ncvreg(x = x, 
                      y = y, 
                      family = "gaussian",
                      standardize = FALSE, 
                      intercept = FALSE)
  # test$beta
  
  r_spsp1 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = MCP.ncvreg,
                        init = 1, standardize = FALSE, intercept = FALSE)
  
  est_beta_1 <- r_spsp1$beta_SPSP
  # colnames(X)[which(est_beta_1 != 0)]
  # head(est_beta_1)
  
  r_spsp5 <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = MCP.ncvreg,
                        init = 5, standardize = FALSE, intercept = FALSE)
  est_beta_5 <- r_spsp5$beta_SPSP
  # head(est_beta_5)
  
  r_spsp5_std <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = MCP.ncvreg,
                            init = 5, standardize = TRUE, intercept = FALSE)
  est_beta_5_std <- r_spsp5_std$beta_SPSP
  # head(est_beta_5_std)
  
  r_spsp5_std_int <- SPSP::SPSP(x = x, y = y, family = "gaussian", fitfun.SP = MCP.ncvreg,
                                init = 5, standardize = TRUE, intercept = TRUE)
  est_beta_5_std_int <- r_spsp5_std_int$beta_SPSP
  # head(est_beta_5_std_int)
  # r_spsp5_std_int$intercept
  
  # Check if coefficients are same as ols
  expect_equal(coef(test_lm_3), est_beta_1[1:3])
  expect_equal(coef(test_lm_3), est_beta_5[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std[1:3])
  expect_equal(coef(test_lm_std_3), est_beta_5_std_int[1:3])
  
  # Check if intercept is correct
  expect_true(is.na(r_spsp1$intercept))
  expect_true(is.na(r_spsp5$intercept))
  expect_true(is.na(r_spsp5_std$intercept))
  expect_true(!is.na(r_spsp5_std_int$intercept))
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

test_that("Check penalty.factor argument in glmnet() for SPSP()", {
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
  x <- as.matrix(Boston2[,-14]); medv <- Boston2[,14] 
  summary(x)
  x[,"crim"] <- log(x[,"crim"])
  x[,"tax"] <- log(x[,"tax"])
  x[,"lstat"] <- log(x[,"lstat"])
  x[,"dis"] <- log(x[,"dis"])
  # x[,"rad"] <- log(x[,"rad"])
  medv <- scale(medv)
  x[,"age"] <- log(x[,"age"])
  
  for (i in 1:(ncol(x))){
    x[,i] <- scale(x[,i])
  }
  # summary(x)
  
  # round(cor(x), 2)
  # round(cor(Boston[,-14]), 2)
  # summary(y)
  
  expect_error(SPSP::SPSP(x = x, y = medv, family = "gaussian", fitfun.SP = lasso.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  test1 <- SPSP::SPSP(x = x, y = medv, family = "gaussian", fitfun.SP = lasso.glmnet,
                      args.fitfun.SP = list(penalty.factor = c(rep(1, 5), 0, rep(1, 6), 0)),
                      init = 1, standardize = T, intercept = FALSE)
  # test1$nonzero

  fit1 <- attr(test1, "glmnet.fit")
  # fit1$beta
  
  # library(plotmo) 
  
  # plot_glmnet(fit1, label=15, xvar ="lambda")
  
  expect_error(SPSP::SPSP(x = x, y = medv, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                          init = 5, standardize = T, intercept = FALSE),
               NA)
  
  test2 <- SPSP::SPSP(x = x, y = medv, family = "gaussian", fitfun.SP = adalassoCV.glmnet,
                      args.fitfun.SP = list(penalty.factor = c(rep(1, 5), 0, rep(1, 6), 0)),
                      init = 1, standardize = T, intercept = FALSE)
  # test2$nonzero
  fit2 <- attr(test2, "glmnet.fit")
  
  # plot_glmnet(fit2, label=15, xvar ="lambda")
  
})
