# Given the same inputs, the solver should produce identical output.

test_that("cpmnet with identical inputs yields identical beta", {
  dat <- sim_cpm(n = 60, p = 4, seed = 700)
  fit1 <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 8,
                 verbose = FALSE)
  fit2 <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 8,
                 verbose = FALSE)
  expect_equal(fit1$beta, fit2$beta)
  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$logL, fit2$logL)
})

test_that("cv.cpmnet sequential with set.seed is deterministic", {
  dat <- sim_cpm(n = 60, p = 4, seed = 701)
  set.seed(1234)
  cv1 <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 6,
                   nfolds = 3, parallel = FALSE, family = "probit")
  set.seed(1234)
  cv2 <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 6,
                   nfolds = 3, parallel = FALSE, family = "probit")
  expect_equal(cv1$lambda.min, cv2$lambda.min)
  expect_equal(cv1$cvm, cv2$cvm)
})

test_that("predict is a pure function of (fit, newx, s, type)", {
  fit <- .small_fit()
  p <- nrow(fit$beta)
  newx <- matrix(rnorm(10 * p), 10, p)
  p1 <- predict(fit, newx = newx, type = "median", s = fit$lambda[5])
  p2 <- predict(fit, newx = newx, type = "median", s = fit$lambda[5])
  expect_equal(p1, p2)
})
