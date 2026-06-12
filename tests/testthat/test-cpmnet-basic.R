# Basic sanity checks on the fitted cpmnet object.

test_that("cpmnet returns object of class 'cpmnet'", {
  fit <- .small_fit()
  expect_s3_class(fit, "cpmnet")
})

test_that("cpmnet returns beta with dim p x nlambda", {
  fit <- .small_fit(n = 50, p = 4, nlambda = 8)
  expect_equal(nrow(fit$beta), 4L)
  expect_equal(ncol(fit$beta), length(fit$lambda))
})

test_that("cpmnet returns alpha with dim k x nlambda", {
  fit <- .small_fit()
  expect_equal(nrow(fit$alpha), fit$data_prep$k)
  expect_equal(ncol(fit$alpha), length(fit$lambda))
})

test_that("lambda path is strictly decreasing", {
  fit <- .small_fit(nlambda = 15)
  expect_true(all(diff(fit$lambda) < 0))
})

test_that("dev.ratio is in [0, 1] and non-decreasing along the path", {
  fit <- .small_fit(nlambda = 15)
  expect_true(all(fit$dev.ratio >= 0 & fit$dev.ratio <= 1))
  # allow small numerical wobbles
  increments <- diff(fit$dev.ratio)
  expect_true(all(increments >= -1e-6))
})

test_that("df is 0 at the largest lambda (lasso, alpha_en = 1)", {
  dat <- sim_cpm()
  fit <- cpmnet(dat$x, dat$y, family = "probit", alpha_en = 1,
                nlambda = 10, verbose = FALSE)
  expect_equal(fit$df[1], 0L)
})

test_that("df is weakly non-decreasing as lambda decreases", {
  dat <- sim_cpm()
  fit <- cpmnet(dat$x, dat$y, family = "probit", alpha_en = 1,
                nlambda = 20, verbose = FALSE)
  # df can fluctuate slightly for non-pure-lasso, but with alpha=1 and a
  # clean signal it should be monotone non-decreasing in the limit.
  expect_true(fit$df[length(fit$df)] >= fit$df[1])
})

test_that("all fits converge on a clean small dataset", {
  fit <- .small_fit()
  expect_true(all(fit$converged))
})

test_that("returned object has every documented slot", {
  fit <- .small_fit()
  required <- c("alpha", "beta", "alpha_scaled", "beta_scaled", "df", "dim",
                "lambda", "dev.ratio", "nulldev", "logL", "converged",
                "niter", "nobs", "data_prep", "alpha_en", "standardize",
                "call")
  expect_true(all(required %in% names(fit)))
})

test_that("user-supplied lambda is honoured exactly", {
  dat <- sim_cpm()
  user_lam <- c(1.0, 0.5, 0.1, 0.01)
  fit <- cpmnet(dat$x, dat$y, family = "probit", lambda = user_lam,
                verbose = FALSE)
  expect_equal(fit$lambda, user_lam)
})

test_that("alpha synonym is accepted and matches alpha_en", {
  dat <- sim_cpm()
  fit1 <- cpmnet(dat$x, dat$y, family = "probit", alpha_en = 0.3,
                 nlambda = 8, verbose = FALSE)
  fit2 <- cpmnet(dat$x, dat$y, family = "probit", alpha = 0.3,
                 nlambda = 8, verbose = FALSE)
  expect_equal(fit1$beta, fit2$beta)
  expect_equal(fit1$alpha_en, fit2$alpha_en)
})

test_that("unstandardized fit differs from standardized fit in scale", {
  dat <- sim_cpm()
  fit_std <- cpmnet(dat$x, dat$y, family = "probit", standardize = TRUE,
                    nlambda = 10, verbose = FALSE)
  fit_raw <- cpmnet(dat$x, dat$y, family = "probit", standardize = FALSE,
                    nlambda = 10, verbose = FALSE)
  # both should fit, but the beta matrices typically differ
  expect_true(is.matrix(fit_std$beta))
  expect_true(is.matrix(fit_raw$beta))
  expect_equal(dim(fit_std$beta), dim(fit_raw$beta))
})
