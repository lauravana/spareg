test_that("Results has right class", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x, y)
  expect_equal(class(spar_res),"spar")
})

test_that("Coef returns vector of correct length", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  expect_equal(length(sparcoef$beta),30)
})

test_that("Coef is more sparse for higher threshold", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  sparcoef2 <- coef(spar_res,
                    nu = spar_res$val_res$nu[which(spar_res$val_res$nu==sparcoef$nu)+1])
  expect_equal(all(which(sparcoef$beta==0) %in% which(sparcoef2$beta==0)),TRUE)
})

example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100, seed = 1234)

test_that("Returned coef and preds are correct for fixed screening and projections", {
  x <- example_data$x
  y <- example_data$y
  xnew <- example_data$xtest

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  spar_res <- spar(x,y,
                   nummods=c(2),
                   inds = list(1:(2*nrow(x)),
                               500+1:(2*nrow(x))),
                   RPMs = list(RP1,RP2))
  sparcoef <- coef(spar_res)
  pred     <- predict(spar_res,xnew=xnew)
  expect_equal(sparcoef$nu,0.002285171,tolerance = 1e-6)
  expect_equal(sparcoef$beta[53],0)
  expect_equal(sparcoef$beta[1],0.125971, tolerance = 1e-6)
  expect_equal(pred[1],20.44922,tolerance = 1e-5)
})

test_that("Returned coef and preds are correct for fixed screening and projections for binomial(logit)", {
  x <- example_data$x
  y <- round(1/(1+exp(-example_data$y)))
  xnew <- example_data$xtest

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  spar_res <- spar(x,y,nummods=c(2),
                   inds = list(1:(2*nrow(x)),500+1:(2*nrow(x))),
                   RPMs = list(RP1,RP2),family = binomial(logit))
  sparcoef <- coef(spar_res)
  pred <- predict(spar_res,xnew=xnew)
  expect_equal(sparcoef$nu,0.009850679 ,tolerance = 1e-6)
  expect_equal(sparcoef$beta[11],0)
  expect_equal(sparcoef$beta[1],0.04795905,tolerance = 1e-6)
  expect_equal(pred[1],0.9749038,tolerance = 1e-5)

})

test_that("Columns with zero sd get ceofficient 0", {
  x <- example_data$x
  x[,c(1,11,111)] <- 2
  y <- example_data$y
  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  expect_equal(sparcoef$beta[c(1,11,111)],c(0,0,0))
})

test_that("Thresholding can be avoided ", {
  x <- example_data$x
  y <- example_data$y
  set.seed(123)
  spar_res <- spar(x, y, screencoef = screen_glmnet(),
                   nus = 0, model = spar_glm())
  sparcoef <- coef(spar_res)
  expect_equal(sparcoef$beta[c(4,10)],c(0,0))
})

test_that("Data splitting delivers different results", {
  x <- example_data$x
  y <- example_data$y
  set.seed(123)
  spar_res <- spar(x,y,
                   screencoef = screen_cor(split_data_prop = 0.25))
  set.seed(123)
  spar_res3 <- spar(x,y,
                    screencoef = screen_cor(split_data = TRUE)) # this should not work
  set.seed(123)
  spar_res2 <- spar(x,y,
                    screencoef = screen_cor())
  sparcoef <- coef(spar_res)
  sparcoef2 <- coef(spar_res2)
  sparcoef3 <- coef(spar_res3)
  expect_true(any(sparcoef$beta[c(8,9,12)] != sparcoef2$beta[c(8,9,12)]))
  expect_equal(sparcoef2$beta[c(8,9,12)], sparcoef3$beta[c(8,9,12)])
})

test_that("Test the gaussian rp", {
  x <- example_data$x
  y <- example_data$y
  set.seed(123)
  spar_g_res <- spar(x,y, screencoef = screen_glmnet(),
                     rp = rp_gaussian())
  set.seed(123)
  spar_g_res2 <- spar(x,y,screencoef = screen_glmnet(),
                      rp = rp_gaussian(sd = 0.1))
  ##
  expect_equal(round(spar_g_res$val_res$Meas[1], 2), 17873.92)
  expect_equal(round(spar_g_res2$val_res$Meas[1], 2), 17873.92)
  expect_equal(spar_g_res$val_res$Meas[1], spar_g_res2$val_res$Meas[1])
})

test_that("Test the sparse rp", {
  x <- example_data$x
  y <- example_data$y
  set.seed(12345)
  spar_sparse_res <- spar(x,y,screencoef = screen_glmnet(),
                          rp = rp_sparse())
  set.seed(12345)
  spar_sparse_res2 <- spar(x,y, screencoef = screen_glmnet(),
                           rp = rp_sparse(psi = 0.01))
  expect_equal(round(spar_sparse_res$val_res$Meas[1], 2), 19028.88)
  expect_equal(round(spar_sparse_res2$val_res$Meas[1], 2), 19004.14)
  expect_true(spar_sparse_res$val_res$Meas[1] != spar_sparse_res2$val_res$Meas[1])
})
test_that("Test the CW rp", {
  x <- example_data$x
  y <- example_data$y
  set.seed(123)
  spar_cw_res <- spar(x,y, screencoef = screen_glmnet(),
                      rp = rp_cw())
  expect_equal(round(spar_cw_res$val_res$Meas[1], 2), 16841.77)
})

test_that("Test the screen_glm() with poisson family", {
  x <- example_data$x
  y <- example_data$y
  yval <- round(abs(y))
  set.seed(123)
  spar_screen_glm <- spar(x,yval, family = poisson(),
                          screencoef = screen_marglik(),
                          rp = rp_gaussian())
  expect_equal(round(spar_screen_glm$val_res$Meas[1], 2), 1431.17)
})

test_that("Get same results with parallel option", {
  x <- example_data$x
  y <- example_data$y
  set.seed(123)
  spar_res <- spar(x, y, screencoef = screen_cor(), rp = rp_gaussian(),
                   seed = 123, set.seed.iteration = TRUE)
  if (requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(2, "PSOCK")
    doParallel::registerDoParallel(cl)
    suppressMessages(
      spar_res2 <- spar(x, y, screencoef = screen_cor(), rp = rp_gaussian(),
                        parallel = TRUE, set.seed.iteration = TRUE, seed = 123))
    parallel::stopCluster(cl)
    expect_equal(spar_res$betas[1:10],  spar_res2$betas[1:10])
  }
})


# Tests expecting errors

test_that("Get errors for data.frame input x", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  expect_error(spar(x,y))
})

test_that("Get errors for categorical input y", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- factor(rep(c("a","b"),5))
  expect_error(spar(x,y))
})

test_that("Get errors for mslow > msup", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  expect_error(spar(x,y,mslow = 15,msup = 10))
})

test_that("Get errors for msup > nscreen", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  expect_error(spar(x,y,nscreen=18, msup = 20))
})

test_that("Get errors for to small length of inds and RPMs lists", {
  x <- example_data$x
  y <- example_data$y

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  expect_error(spar(x,y,nummods=c(3),inds = list(1:(2*nrow(x)),500+1:(2*nrow(x))),
                    RPMs = list(RP1,RP2)))
})

test_that("Get errors for prediction when xnew has wrong dimensions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar(x,y)
  xnew <- example_data$xtest
  expect_error(predict(spar_res,xnew=xnew[,-1]))
})

test_that("Get errors for classification validation measure for non-binomial family", {
  x <- example_data$x
  y <- example_data$y
  expect_error(spar(x,y,measure = "1-auc"))
})
