# Edge cases: extreme n/p ratios, sparse signals, correlated predictors, binary
# and single-valued y, constant columns, weights, penalty factors, a singular
# intercept Hessian and the first fit of a path. Each test is kept small so the
# suite stays fast.

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
  skip_if_not_installed("MASS")
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

dat_rf <- function(n = 80, p = 5, seed = 21) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p)
  list(x = x, y = exp(x[, 1] - 0.5 * x[, 2] + rnorm(n)))
}

test_that("a constant column gets coefficient 0 and leaves the rest of the fit as without it", {
  d <- dat_rf()
  xc <- cbind(d$x[, 1:2], 0, d$x[, 3:5])
  for (core in c("group", "elementwise")) {
    old <- options(cpmnet.shuffle_seed = 3L); on.exit(options(old), add = TRUE)
    f0 <- cpmnet(d$x, d$y, core = core, nlambda = 10, verbose = FALSE)
    f1 <- cpmnet(xc, d$y, core = core, nlambda = 10, verbose = FALSE)
    expect_equal(nrow(f1$beta), 6L)
    expect_true(all(f1$beta[3, ] == 0))
    expect_equal(f1$beta[-3, ], f0$beta, tolerance = 1e-12)
    expect_equal(f1$constant_columns, 3L)
    expect_equal(ncol(predict(f1, newx = xc[1:3, ])), 10L)
  }
  expect_warning(cpmnet(xc, d$y, nlambda = 5), "constant column")
})

test_that("cross-validation survives a fold whose training rows hold a constant column", {
  d <- dat_rf(n = 60)
  rare <- rep(0, 60); rare[c(1, 6)] <- 1
  x <- cbind(d$x, rare)
  for (core in c("group", "elementwise")) {
    cv <- cv.cpmnet(x, d$y, foldid = rep(1:5, length.out = 60), nlambda = 8,
                    parallel = FALSE, core = core, verbose = FALSE)
    expect_equal(cv$nfolds_successful, 5L)
  }
})

test_that("every column constant and a single-valued y stop with a plain message", {
  d <- dat_rf()
  expect_error(cpmnet(matrix(1, 80, 3), d$y), "every column of x is constant")
  expect_error(cpmnet(d$x, rep(2, 80)), "single value")
})

test_that("the two cores fit the same model under unequal weights", {
  d <- dat_rf(n = 150)
  set.seed(4); w <- runif(150, 0.2, 3)
  fe <- cpmnet(d$x, d$y, weights = w, core = "elementwise", nlambda = 15, verbose = FALSE)
  fg <- cpmnet(d$x, d$y, weights = w, core = "group", lambda = fe$lambda,
               tol = 1e-10, verbose = FALSE)
  expect_equal(fg$beta, fe$beta, tolerance = 1e-5)
})

test_that("cross-validation runs with weights, and unit weights give the unweighted result", {
  d <- dat_rf(n = 100)
  fold <- rep(1:5, length.out = 100)
  old <- options(cpmnet.shuffle_seed = 5L); on.exit(options(old), add = TRUE)
  set.seed(1); a <- cv.cpmnet(d$x, d$y, foldid = fold, nlambda = 8, parallel = FALSE)
  set.seed(1); b <- cv.cpmnet(d$x, d$y, foldid = fold, nlambda = 8, parallel = FALSE,
                              weights = rep(1, 100))
  expect_equal(b$cvm, a$cvm, tolerance = 1e-10)
  set.seed(2); w <- runif(100, 0.5, 2)
  cw <- cv.cpmnet(d$x, d$y, foldid = fold, nlambda = 8, parallel = FALSE, weights = w)
  expect_true(all(is.finite(cw$cvm)))
})

test_that("with penalty factors, every penalized coefficient is zero at the first lambda", {
  d <- dat_rf()
  pf <- c(0.2, rep(1, 4))
  for (core in c("group", "elementwise")) {
    f <- cpmnet(d$x, d$y, penalty_factor = pf, core = core, nlambda = 10, verbose = FALSE)
    expect_true(all(f$beta[, 1] == 0))
    expect_true(any(f$beta[, 2] != 0))
  }
})

test_that("a grouping with an unpenalized column accepts one factor per penalized group", {
  set.seed(7); n <- 120
  fac <- factor(sample(letters[1:3], n, TRUE)); xf <- model.matrix(~ fac)[, -1]
  x <- cbind(xf, matrix(rnorm(n * 3), n, 3), z = rnorm(n))
  y <- exp(xf[, 1] + x[, 3] + rnorm(n))
  f <- cpmnet(x, y, group = c(1, 1, 2, 3, 4, 0), penalty_factor = c(2, 1, 1, 1),
              nlambda = 8, verbose = FALSE)
  expect_s3_class(f, "cpmnet")
  expect_true(all(f$beta[6, ] != 0))
})

test_that("the grouped core honours a supplied lambda_max", {
  d <- dat_rf()
  f <- cpmnet(d$x, d$y, lambda_max = 0.7, nlambda = 5, verbose = FALSE)
  expect_equal(f$lambda[1], 0.7)
})

test_that("a singular intercept Hessian keeps the lambdas fitted before it", {
  set.seed(1)
  x <- matrix(rnorm(50 * 5), 50, 5)
  y <- x[, 1] * 10 + rnorm(50, sd = 0.01)
  for (core in c("group", "elementwise")) {
    f <- suppressWarnings(cpmnet(x, y, core = core, nlambda = 20))
    expect_s3_class(f, "cpmnet")
    expect_lt(length(f$lambda), 20L)
    expect_true(all(is.finite(f$beta)))
  }
  expect_warning(cpmnet(x, y, nlambda = 20), "truncated")
})

test_that("the conditional CDF at thresholds equals a direct calculation", {
  d <- dat_rf()
  f <- cpmnet(d$x, d$y, family = "logistic", nlambda = 6, verbose = FALSE)
  ym <- f$data_prep$y_mapping
  cuts <- c(min(ym) - 1, ym[3], (ym[10] + ym[11]) / 2, max(ym), max(ym) + 1)
  got <- cpmnet:::predict_cdf_cpmnet(f, d$x[1:4, ], thresholds = cuts)
  for (l in 1:6) for (t in seq_along(cuts)) {
    j <- findInterval(cuts[t], ym)
    want <- if (j == 0) rep(0, 4) else if (j == length(ym)) rep(1, 4) else
      1 - plogis(f$alpha[j, l] + drop(d$x[1:4, ] %*% f$beta[, l]))
    expect_equal(got[, l, t], want, tolerance = 1e-12)
  }
})

test_that("the first fit of a computed path has every coefficient exactly zero", {
  set.seed(9); n <- 120
  x <- matrix(rnorm(n * 8), n, 8); y <- exp(0.3 * x[, 1] - 0.2 * x[, 2] + rnorm(n))
  for (core in c("group", "elementwise")) for (fam in c("logistic", "probit", "loglog", "cloglog", "cauchit"))
    for (a in c(1, 0.5, 0.3)) {
      f <- cpmnet(x, y, family = fam, alpha_en = a, core = core, nlambda = 5, verbose = FALSE)
      expect_true(all(f$beta[, 1] == 0), info = paste(core, fam, a))
      expect_equal(f$df[1], 0)
    }
  # ridge has no lambda at which every coefficient is zero, so its first fit is left alone
  f <- cpmnet(x, y, alpha_en = 0, nlambda = 5, verbose = FALSE)
  expect_true(all(f$beta[, 1] != 0))
})
