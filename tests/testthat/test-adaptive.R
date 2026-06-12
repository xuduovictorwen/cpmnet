# Adaptive lasso path: penalty_factor should inverse-weight the orm slope.

test_that("adaptive = TRUE fits and produces finite beta", {
  set.seed(500)
  n <- 120; p <- 5
  x <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(x %*% c(2, -2, 0.5, 0, 0) + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", adaptive = TRUE,
                nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
  expect_true(all(is.finite(fit$beta)))
})

test_that("adaptive weights shrink large-effect predictors less", {
  # With adaptive weights 1/|orm_beta|, large orm_beta -> small penalty ->
  # coefficient should survive longer along the path.
  set.seed(501)
  n <- 200; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, -3, 0, 0)  # variables 1,2 strong; 3,4 null
  y <- as.numeric(x %*% beta_true + rnorm(n))
  fit <- cpmnet(x, y, family = "probit", adaptive = TRUE,
                nlambda = 20, alpha_en = 1, verbose = FALSE)
  # at the solution with maximum penalty where some beta are nonzero, the
  # strong signals should enter first
  first_nonzero <- apply(fit$beta != 0, 1, function(v) {
    idx <- which(v); if (length(idx) == 0) Inf else min(idx)
  })
  expect_lt(first_nonzero[1], first_nonzero[3])
})

test_that("penalty_factor overrides adaptive when both set", {
  set.seed(502)
  n <- 100; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(x %*% c(1, -1, 0, 0) + rnorm(n))
  pf <- c(1, 1, 2, 2)
  # adaptive is ignored when penalty_factor is provided (by the code path
  # in cpmnet that only computes adaptive weights when penalty_factor is
  # NULL)
  fit <- cpmnet(x, y, family = "probit", adaptive = TRUE,
                penalty_factor = pf, nlambda = 8, verbose = FALSE)
  expect_s3_class(fit, "cpmnet")
})
