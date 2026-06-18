example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100, seed = 1234)

test_that("Results has right class", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x, y, rp = rp_gaussian(), nfolds = 4L,
                      model = spar_glm())
  expect_equal(class(spar_res),"spar.cv")
})

test_that("Coef returns vector of correct length", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x, y, model = spar_glm(),
                      nfolds = 4L)
  sparcoef <- coef(spar_res)
  expect_equal(length(sparcoef$beta),30)
})

test_that("No screening: fast_fit = fix_rpm delivers different results than fix_rpm_and_inds but same threshold", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x, y, rp = rp_gaussian(),
                      model = spar_glm(),
                      nfolds = 4L, seed = 1234)
  spar_res2 <- spar.cv(x, y, rp = rp_gaussian(),
                       model = spar_glm(),
                       nfolds = 4L, seed = 1234,
                       fast_fit = "fix_rpm")
  a <- get_model(spar_res, "best",, x = x, y = y)
  b <- get_model(spar_res2, "best", x = x, y = y)
  expect_true(get_measure(a)[,3] == get_measure(b)[,3])
  expect_true(get_measure(a)[,1] == get_measure(b)[,1])
  expect_error(coef(spar_res2))
  expect_error(predict(spar_res2, x))
})


test_that("With screening: fast_fit = fix_rpm delivers different results than fix_rpm_and_inds but same threshold", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x, y, rp = rp_gaussian(),
                      screencoef = screen_cor(),
                      model = spar_glm(),
                      nfolds = 4L, seed = 1234)
  spar_res2 <- spar.cv(x, y, rp = rp_gaussian(),
                       model = spar_glm(),
                       screencoef = screen_cor(),
                       nfolds = 4L, seed = 1234,
                       fast_fit = "fix_rpm")

  a <- get_model(spar_res, "best")
  b <- get_model(spar_res2, "best", x = x, y = y)
  expect_false(get_measure(a)[,3] == get_measure(b)[,3])
  expect_error(coef(spar_res2))
  expect_error(predict(spar_res2, x))
})

test_that("No screening: fast_fit = fix_rpm delivers same results as none due to seeds ", {
  x <- matrix(rnorm(50*60), ncol = 60)
  y <- rnorm(50)
  spar_res <- spar.cv(x, y, rp = rp_gaussian(),
                      model = spar_glm(),
                      nfolds = 4L, seed = 1234,
                      fast_fit = "none")

  spar_res2 <- spar.cv(x, y, rp = rp_gaussian(),
                       model = spar_glm(),
                       nfolds = 4L, seed = 1234,
                       fast_fit = "fix_rpm")
  a <- get_model(spar_res, "best", x = x, y = y)
  b <- get_model(spar_res2, "best", x = x, y = y)
  expect_true(get_measure(a)[,3] == get_measure(b)[,3])
  expect_error(coef(spar_res))
  expect_error(predict(spar_res, x))
})

test_that("No screening: fast_fit = fix_rpm delivers diff results than none ", {
  x <- matrix(rnorm(50*60), ncol = 60)
  xnew <- matrix(rnorm(50*60), ncol = 60)
  y <- rnorm(50)
  ynew <- rnorm(50)
  nu <- exp(seq(-3, -1, length.out = 20))
  spar_res <- spar.cv(x, y, rp = rp_gaussian(),
                      model = spar_glm(),
                      nus = nu,
                      screencoef = screen_cor(nscreen = 30),
                      nfolds = 4L, seed = 1234,
                      fast_fit = "none")
  spar_res2 <- spar.cv(x, y, rp = rp_gaussian(),
                       model = spar_glm(),
                       nus = nu,
                       screencoef = screen_cor(nscreen = 30),
                       nfolds = 4L, seed = 1234,
                       fast_fit = "fix_rpm")
  spar_res3 <- spar.cv(x, y, rp = rp_gaussian(),
                       model = spar_glm(),
                       nus = nu,
                       screencoef = screen_cor(nscreen = 30),
                       nfolds = 4L, seed = 1234)
  a <- get_model(spar_res, "best", x = x, y = y,
                 xval = xnew,
                 yval = ynew)
  b <- get_model(spar_res2, "best", x = x, y = y,
                 xval = xnew,
                 yval = ynew)
  cc <- get_model(spar_res3, "best", x = x, y = y,
                  xval = xnew,
                  yval = ynew)
  expect_true(get_measure(a)[, 3] <= get_measure(cc)[,3])
})


test_that("Test get_intercept(), get_coef(), get_model() extractors and that coef is more sparse for 1se rule", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar.cv(x, y, rp = rp_gaussian(),  nfolds = 2L, seed = 123)
  sparcoef <- coef(spar_res, opt_par = "best")
  sparcoef2 <- coef(spar_res,opt_par = "1se")
  a <- get_model(spar_res, "best")
  b <- get_model(spar_res, "1se")
  expect_equal(unname(get_intercept(sparcoef)), 1.824784, tolerance = 1e-5)
  expect_equal(unname(get_coef(sparcoef)[1]), 0.2364498, tolerance = 1e-5)
  expect_equal(all(which(sparcoef$beta==0) %in% which(sparcoef2$beta==0)),TRUE)
  expect_equal(nrow(a$val_res), 1)
  expect_equal(nrow(b$val_res), 1)
})

test_that("Validated nu and nummod values are same as the ones for initial SPAR fit", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,screencoef = screen_glmnet(),
                      nummods=c(10,15), model = spar_glm(),
                      nfolds = 4L)
  expect_equal(unique(spar_res$val_res$nu),as.numeric(spar_res$nus))
  expect_equal(unique(spar_res$val_res$nummod),as.numeric(spar_res$nummods))

})

test_that("Columns with zero sd get coefficient 0", {
  x <- example_data$x
  x[,c(1,11,111)] <- 2
  y <- example_data$y
  spar_res <- spar.cv(x, y, screencoef = screen_cor(),
                      measure = "mae", model = spar_glm(),
                      nfolds = 4L)
  sparcoef <- coef(spar_res, opt_par = "best")
  sparcoef2 <- coef(spar_res, opt_par = "1se")
  expect_equal(unname(sparcoef$beta[c(1,11,111)]),c(0,0,0))
  expect_equal(unname(sparcoef2$beta[c(1,11,111)]),c(0,0,0))

})

test_that("Columns with zero sd get coefficient 0 when fast_fit = none", {
  x <- example_data$x
  x[,c(1,11,111)] <- 2
  y <- example_data$y
  spar_res <- spar.cv(x, y, screencoef = screen_cor(),
                      measure = "mae", model = spar_glm(),
                      nfolds = 4L, fast_fit = "none")
  a <- get_model(spar_res, "best", x = x, y = y)
  sparcoef <- coef(a, opt_par = "best")
  sparcoef2 <- coef(a)
  expect_equal(unname(sparcoef$beta[c(1,11,111)]),c(0,0,0))
  expect_equal(unname(sparcoef2$beta[c(1,11,111)]),c(0,0,0))
})

test_that("Test get_measure() extractor and predictions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar.cv(x, y, nfolds = 2L, seed = 123)
  a <- get_model(spar_res, "best")
  b <- get_model(spar_res, "1se")
  expect_equal(nrow(get_measure(a)), 1)
  expect_equal(nrow(get_measure(b)), 1)
  expect_equal(get_measure(a)$deviance, 276.0395, tolerance = 1e-5)
  expect_equal(get_measure(a)$numactive, 311, tolerance = 1e-5)
  expect_equal(get_measure(b)$deviance, 920.7225, tolerance = 1e-5)
  expect_equal(get_measure(b)$numactive, 252, tolerance = 1e-5)
  xnew <- example_data$xtest
  pred  <- predict(spar_res,xnew = xnew)
  pred2 <- predict(spar_res,xnew = xnew, opt_par = "best")
  pred3 <- predict(spar_res, xnew = xnew, opt_par = "1se")
  pred4 <- predict(spar_res, xnew = xnew, opt_par = "1se", avg_type = "link")
  expect_equal(length(pred), nrow(xnew))
  expect_equal(pred[1], pred2[1], tolerance = 1e-5)
  expect_equal(pred3[1], pred4[1], tolerance = 1e-5)
  expect_lt(pred3[1], pred2[1])
})

# Tests expecting errors

test_that("Get errors for input x not data.frame or matrix", {
  x <- list("1"=1:10,"2"=(-1)^(1:12),"3"=rnorm(12),
            "4"=rnorm(12),"5"=runif(12),"6"=runif(12))
  y <- rnorm(6)
  expect_error(spar.cv(x,y,nfolds=2L, model = spar_glm()))
})

test_that("Get errors for input x non-numeric data.frame", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  x$X1[1] <- "a"
  y <- rnorm(10)
  expect_error(spar.cv(x,y))
})

test_that("Get errors for categorical input y", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- factor(rep(c("a","b"),5))
  expect_error(spar.cv(x,y))
})

test_that("Get errors for prediction when xnew has wrong dimensions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar.cv(x,y,model = spar_glm(),
                      nfolds = 3L)
  xnew <- example_data$xtest
  expect_error(predict(spar_res,xnew=xnew[,-1]))
})

test_that("Get errors for classification validation measure for non-binomial family", {
  x <- example_data$x
  y <- example_data$y
  expect_error(spar.cv(x,y,measure = "1-auc",model = spar_glm()))
})




