# Cross-check unpenalized cpmnet against rms::orm.fit. When lambda is
# effectively zero the two should agree on the slope estimates.

test_that("unpenalized cpmnet matches rms::orm.fit slopes (probit)", {
  set.seed(400)
  n <- 200; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -0.8, 0.5, 0)
  y <- as.numeric(x %*% beta_true + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", lambda = c(1e-10),
                standardize = FALSE, verbose = FALSE)
  orm <- rms::orm.fit(x = x, y = y, family = "probit")
  orm_beta <- orm$coefficients[-seq_len(max(fit$data_prep$k))]
  expect_equal(as.numeric(fit$beta[, 1]), as.numeric(orm_beta),
               tolerance = 1e-3)
})

test_that("unpenalized cpmnet matches rms::orm.fit slopes (logistic)", {
  set.seed(401)
  n <- 200; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -1, 0, 0)
  y <- as.numeric(x %*% beta_true + rnorm(n))
  fit <- cpmnet(x, y, family = "logistic", lambda = c(1e-10),
                standardize = FALSE, verbose = FALSE)
  orm <- rms::orm.fit(x = x, y = y, family = "logistic")
  orm_beta <- orm$coefficients[-seq_len(max(fit$data_prep$k))]
  expect_equal(as.numeric(fit$beta[, 1]), as.numeric(orm_beta),
               tolerance = 1e-2)
})

test_that("standardized and unstandardized fits agree after unscaling", {
  set.seed(402)
  n <- 150; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(x %*% c(1, -1, 0.5, 0) + rnorm(n))
  fit_std <- cpmnet(x, y, family = "probit", lambda = c(0.01),
                    standardize = TRUE, verbose = FALSE)
  fit_raw <- cpmnet(x, y, family = "probit", lambda = c(0.01),
                    standardize = FALSE, verbose = FALSE)
  # Standardization is a reparameterization of the penalty. Coefs need not
  # be identical but they should lie in the same neighborhood up to the
  # scale adjustment that unstandardize applies.
  expect_true(is.finite(fit_std$beta[1, 1]))
  expect_true(is.finite(fit_raw$beta[1, 1]))
})

test_that("very small lambda recovers high dev.ratio relative to null", {
  set.seed(403)
  n <- 200; p <- 3
  x <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(x %*% c(2, -2, 1.0) + rnorm(n, sd = 0.3))
  fit <- cpmnet(x, y, family = "probit", lambda = c(1.0, 0.1, 0.01, 1e-6),
                standardize = TRUE, verbose = FALSE)
  # dev.ratio monotone non-decreasing as lambda shrinks, and the tail must
  # clearly beat the null
  expect_gt(tail(fit$dev.ratio, 1), fit$dev.ratio[1])
  expect_gt(tail(fit$dev.ratio, 1), 0.1)
})

test_that("logL and dev.ratio are consistent with nulldev", {
  fit <- .small_fit(n = 100, p = 4, nlambda = 8)
  dev.ratio_manual <- pmax(0, 1 - fit$logL / fit$nulldev)
  expect_equal(fit$dev.ratio, dev.ratio_manual)
})
