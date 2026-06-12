# Edge cases: extreme n/p ratios, sparse signals, correlated predictors,
# binary y. Each test is kept small so the suite stays fast.

test_that("n > p fits cleanly", {
  dat <- sim_cpm(n = 200, p = 5)
  fit <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 10,
                verbose = FALSE)
  expect_true(all(fit$converged))
})

test_that("n < p fits without error", {
  set.seed(7)
  n <- 30; p <- 50
  x <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -1, 0.5, rep(0, p - 3))
  y <- as.numeric(x %*% beta_true + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", nlambda = 10, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_equal(nrow(fit$beta), p)
})

test_that("n == p fits", {
  set.seed(8)
  n <- 25; p <- 25
  x <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(x %*% c(1, rep(0, p - 1)) + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", nlambda = 10, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
})

test_that("binary y is accepted and produces non-NA output", {
  set.seed(9)
  x <- matrix(rnorm(60 * 4), 60, 4)
  y <- rbinom(60, 1, 0.5)
  fit <- cpmnet(x, y, family = "logistic", nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_true(all(is.finite(fit$beta)))
})

test_that("heavily tied y is handled via rank mapping", {
  set.seed(10)
  x <- matrix(rnorm(80 * 4), 80, 4)
  y <- sample(1:5, 80, replace = TRUE)  # many ties
  fit <- cpmnet(x, y, family = "probit", nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_equal(fit$data_prep$k, 4L)  # 5 unique levels - 1
})

test_that("correlated predictors (rho = 0.9) fit without error", {
  set.seed(11)
  n <- 80; p <- 5
  Sigma <- 0.9^abs(outer(1:p, 1:p, "-"))
  x <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  y <- as.numeric(x %*% c(1, -1, rep(0, p - 2)) + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", nlambda = 10, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_true(all(fit$converged))
})

test_that("near-constant predictor with tiny variance still fits", {
  set.seed(12)
  x <- matrix(rnorm(60 * 4), 60, 4)
  x[, 4] <- x[, 4] * 1e-6  # nearly constant after standardization
  y <- as.numeric(x[, 1] + rnorm(60))
  fit <- cpmnet(x, y, family = "probit", nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
})

test_that("single predictor (p = 1) fits", {
  set.seed(13)
  n <- 60
  x <- matrix(rnorm(n), n, 1)
  y <- as.numeric(x + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_equal(nrow(fit$beta), 1L)
})

test_that("constant y triggers a graceful error or warning", {
  set.seed(14)
  x <- matrix(rnorm(40 * 3), 40, 3)
  y <- rep(5, 40)
  expect_error(suppressWarnings(cpmnet(x, y, family = "probit", nlambda = 5,
                                       verbose = FALSE)))
})

test_that("zero-weight rows are dropped transparently", {
  dat <- sim_cpm(n = 80, p = 4)
  w <- rep(1.0, 80); w[1:10] <- 0
  fit <- cpmnet(dat$x, dat$y, family = "probit", weights = w,
                nlambda = 8, verbose = FALSE)
  expect_equal(fit$nobs, 70L)
})
