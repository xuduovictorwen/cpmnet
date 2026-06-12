# Internal predict_* helpers that feed into cv.cpmnet. These are not
# exported, so we access them via the package namespace.

.pq <- cpmnet:::predict_quantile_cpmnet
.pc <- cpmnet:::predict_cdf_cpmnet
.pp <- cpmnet:::predict_psr_cpmnet

test_that("predict_quantile_cpmnet returns n x nlambda x n_tau array", {
  fit <- .small_fit(nlambda = 8, n = 60, p = 4)
  newx <- matrix(rnorm(5 * 4), 5, 4)
  tau <- c(0.25, 0.5, 0.75)
  out <- .pq(fit, newx = newx, s = NULL, tau_levels = tau)
  expect_equal(dim(out), c(5L, length(fit$lambda), 3L))
})

test_that("predict_quantile at tau = 0.5 approximately matches median", {
  fit <- .small_fit(nlambda = 8, n = 60, p = 4)
  newx <- matrix(rnorm(5 * 4), 5, 4)
  qpred <- .pq(fit, newx = newx, s = fit$lambda[4], tau_levels = 0.5)
  mpred <- predict(fit, newx = newx, type = "median", s = fit$lambda[4])
  # both walk the CDF step function; equality should hold exactly on the
  # ordered support
  expect_equal(as.numeric(qpred[, 1, 1]), as.numeric(mpred), tolerance = 1e-8)
})

test_that("predict_quantile rejects tau outside (0, 1)", {
  fit <- .small_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_error(.pq(fit, newx = newx, tau_levels = c(0.1, 1.0)),
               "tau_levels must lie strictly in \\(0, 1\\)")
})

test_that("predict_quantile rejects unsorted tau", {
  fit <- .small_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_error(.pq(fit, newx = newx, tau_levels = c(0.9, 0.1)),
               "tau_levels must be sorted in increasing order")
})

test_that("predict_cdf returns array n x nlambda x n_thresh", {
  fit <- .small_fit(nlambda = 8, n = 60, p = 4)
  newx <- matrix(rnorm(5 * 4), 5, 4)
  th <- quantile(sim_cpm(seed = 1)$y, probs = c(0.25, 0.5, 0.75))
  out <- .pc(fit, newx = newx, s = NULL, thresholds = th)
  expect_equal(dim(out), c(5L, length(fit$lambda), 3L))
})

test_that("predict_cdf values lie in [0, 1]", {
  fit <- .small_fit(nlambda = 8, n = 60, p = 4)
  newx <- matrix(rnorm(5 * 4), 5, 4)
  th <- c(-1, 0, 1)
  out <- .pc(fit, newx = newx, s = NULL, thresholds = th)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("predict_cdf is monotone non-decreasing in threshold per lambda", {
  fit <- .small_fit(nlambda = 8, n = 60, p = 4)
  newx <- matrix(rnorm(5 * 4), 5, 4)
  th <- sort(c(-2, -1, 0, 1, 2))
  out <- .pc(fit, newx = newx, s = NULL, thresholds = th)
  # for each obs and each lambda, cdf should be monotone in threshold
  diffs <- apply(out, c(1, 2), function(v) min(diff(v)))
  expect_true(all(diffs >= -1e-8))
})

test_that("predict_psr_cpmnet returns n x nlambda matrix in [0, 1]", {
  dat <- sim_cpm(n = 60, p = 4, seed = 300)
  fit <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 8,
                verbose = FALSE)
  newx <- dat$x[1:5, , drop = FALSE]
  y_test <- dat$y[1:5]
  out <- .pp(fit, newx = newx, y_test = y_test, s = NULL)
  expect_equal(dim(out), c(5L, length(fit$lambda)))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("predict_cdf rejects non-numeric thresholds", {
  fit <- .small_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_error(.pc(fit, newx = newx, s = NULL, thresholds = "foo"),
               "thresholds must be a non-empty numeric vector")
})

test_that("predict_psr_cpmnet rejects mismatched y_test length", {
  dat <- sim_cpm(n = 40, p = 4, seed = 301)
  fit <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 6,
                verbose = FALSE)
  newx <- dat$x[1:5, , drop = FALSE]
  expect_error(.pp(fit, newx = newx, y_test = dat$y[1:3], s = NULL),
               "y_test must be a numeric vector of length nrow\\(newx\\)")
})
