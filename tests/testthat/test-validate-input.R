# Validation of user input. Most tests drive validation through the public
# cpmnet() / cv.cpmnet() / predict.cpmnet() entry points so that we cover the
# actual error paths users will hit.

test_that("x with NA triggers finite-value error", {
  dat <- sim_cpm()
  dat$x[1, 1] <- NA
  expect_error(cpmnet(dat$x, dat$y, family = "probit", nlambda = 5,
                      verbose = FALSE),
               "x must contain no NA/NaN/Inf")
})

test_that("x with Inf triggers finite-value error", {
  dat <- sim_cpm()
  dat$x[2, 3] <- Inf
  expect_error(cpmnet(dat$x, dat$y, family = "probit", nlambda = 5,
                      verbose = FALSE),
               "x must contain no NA/NaN/Inf")
})

test_that("character x fails numeric check", {
  x <- matrix(letters[1:20], 4, 5)
  y <- 1:4
  expect_error(cpmnet(x, y, family = "probit", nlambda = 5, verbose = FALSE),
               "x must be a numeric matrix")
})

test_that("y with NA triggers finite-value error", {
  dat <- sim_cpm()
  dat$y[1] <- NA
  expect_error(cpmnet(dat$x, dat$y, family = "probit", nlambda = 5,
                      verbose = FALSE),
               "y must contain no NA/NaN/Inf")
})

test_that("non-numeric y fails", {
  dat <- sim_cpm()
  y_chr <- as.character(dat$y)
  expect_error(cpmnet(dat$x, y_chr, family = "probit", nlambda = 5,
                      verbose = FALSE),
               "y must be a numeric vector")
})

test_that("length(y) != nrow(x) fails with sizes", {
  dat <- sim_cpm(n = 20)
  expect_error(cpmnet(dat$x, dat$y[1:10], family = "probit", nlambda = 5,
                      verbose = FALSE),
               "length\\(y\\) = 10 must equal nrow\\(x\\) = 20")
})

test_that("empty x fails", {
  x <- matrix(numeric(0), 0, 0)
  y <- numeric(0)
  expect_error(cpmnet(x, y, family = "probit", nlambda = 5, verbose = FALSE),
               "x must be non-empty")
})

test_that("alpha_en out of range [0, 1] fails with value", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", alpha_en = -0.1,
                      nlambda = 5, verbose = FALSE),
               "alpha_en must be in \\[0, 1\\]; got -0.1")
  expect_error(cpmnet(dat$x, dat$y, family = "probit", alpha_en = 1.5,
                      nlambda = 5, verbose = FALSE),
               "alpha_en must be in \\[0, 1\\]; got 1.5")
})

test_that("non-scalar alpha_en fails", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", alpha_en = c(0.1, 0.2),
                      nlambda = 5, verbose = FALSE),
               "alpha_en must be a single finite numeric value")
})

test_that("unknown family fails and lists valid families", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "banana", nlambda = 5,
                      verbose = FALSE),
               "family must be one of:.*got 'banana'")
})

test_that("family as vector fails", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = c("probit", "logistic"),
                      nlambda = 5, verbose = FALSE),
               "family must be a single character string")
})

test_that("negative weights rejected", {
  dat <- sim_cpm()
  w <- rep(1.0, length(dat$y)); w[3] <- -0.5
  expect_error(cpmnet(dat$x, dat$y, family = "probit", weights = w,
                      nlambda = 5, verbose = FALSE),
               "weights must be non-negative")
})

test_that("all-zero weights rejected", {
  dat <- sim_cpm()
  w <- rep(0.0, length(dat$y))
  expect_error(cpmnet(dat$x, dat$y, family = "probit", weights = w,
                      nlambda = 5, verbose = FALSE),
               "weights must not all be zero")
})

test_that("wrong-length weights rejected", {
  dat <- sim_cpm()
  w <- rep(1.0, length(dat$y) - 1)
  expect_error(cpmnet(dat$x, dat$y, family = "probit", weights = w,
                      nlambda = 5, verbose = FALSE),
               "length\\(weights\\) = .* must equal length\\(y\\)")
})

test_that("non-decreasing lambda sequence warns and sorts", {
  dat <- sim_cpm()
  expect_warning(
    fit <- cpmnet(dat$x, dat$y, family = "probit",
                  lambda = c(0.01, 0.05, 0.1), verbose = FALSE),
    "lambda is not in decreasing order"
  )
  expect_true(all(diff(fit$lambda) < 0))
})

test_that("negative lambda rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      lambda = c(0.1, -0.01), verbose = FALSE),
               "lambda must be non-negative")
})

test_that("nlambda < 2 rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", nlambda = 1,
                      verbose = FALSE),
               "nlambda must be a single integer >= 2")
})

test_that("non-integer nlambda rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", nlambda = 5.5,
                      verbose = FALSE),
               "nlambda must be a single integer >= 2")
})

test_that("max_iter <= 0 rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", max_iter = 0,
                      nlambda = 5, verbose = FALSE),
               "max_iter must be a single positive integer")
})

test_that("tol <= 0 rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit", tol = 0,
                      nlambda = 5, verbose = FALSE),
               "tol must be a single positive number")
})

test_that("lambda_min_ratio < 0 rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      lambda_min_ratio = -0.01, nlambda = 5, verbose = FALSE),
               "lambda_min_ratio must be a single non-negative number")
})

test_that("penalty_factor length mismatch rejected", {
  dat <- sim_cpm()
  bad_pf <- rep(1.0, dat$p - 1)
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      penalty_factor = bad_pf, nlambda = 5, verbose = FALSE),
               "length\\(penalty_factor\\) = .* must equal p = ")
})

test_that("negative penalty_factor rejected", {
  dat <- sim_cpm()
  pf <- rep(1.0, dat$p); pf[2] <- -0.5
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      penalty_factor = pf, nlambda = 5, verbose = FALSE),
               "penalty_factor must be non-negative")
})

test_that("all-zero penalty_factor rejected", {
  dat <- sim_cpm()
  pf <- rep(0.0, dat$p)
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      penalty_factor = pf, nlambda = 5, verbose = FALSE),
               "penalty_factor must not be all zero")
})

test_that("non-logical standardize rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      standardize = "yes", nlambda = 5, verbose = FALSE),
               "standardize must be a single TRUE/FALSE")
})

test_that("non-logical verbose rejected", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      nlambda = 5, verbose = NA),
               "verbose must be a single TRUE/FALSE")
})

test_that("unused dots argument raises", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      nlambda = 5, verbose = FALSE, foo = 3),
               "Unused arguments: foo")
})

test_that("alpha and alpha_en both supplied raises", {
  dat <- sim_cpm()
  expect_error(cpmnet(dat$x, dat$y, family = "probit",
                      alpha_en = 0.5, alpha = 0.5, nlambda = 5,
                      verbose = FALSE),
               "Both `alpha` and `alpha_en` supplied")
})

test_that("cv.cpmnet rejects nfolds < 2", {
  dat <- sim_cpm()
  expect_error(cv.cpmnet(dat$x, dat$y, nfolds = 1, nlambda = 5,
                         parallel = FALSE),
               "nfolds must be at least 2")
})

test_that("cv.cpmnet rejects nfolds > n", {
  dat <- sim_cpm(n = 20)
  expect_error(cv.cpmnet(dat$x, dat$y, nfolds = 25, nlambda = 5,
                         parallel = FALSE),
               "nfolds = 25 exceeds n = 20")
})

test_that("cv.cpmnet rejects malformed foldid (gap)", {
  dat <- sim_cpm(n = 20)
  fid <- rep(c(1L, 3L), 10)  # missing 2
  expect_error(cv.cpmnet(dat$x, dat$y, foldid = fid, nfolds = 2,
                         nlambda = 5, parallel = FALSE),
               "foldid values must form the set 1:max\\(foldid\\)")
})

test_that("cv.cpmnet rejects wrong-length foldid", {
  dat <- sim_cpm(n = 20)
  expect_error(cv.cpmnet(dat$x, dat$y, foldid = rep(1L, 5),
                         nfolds = 2, nlambda = 5, parallel = FALSE),
               "foldid must be an integer vector of length n = 20")
})

test_that("cv.cpmnet rejects tau_levels outside (0, 1)", {
  dat <- sim_cpm()
  expect_error(cv.cpmnet(dat$x, dat$y, type = "pinball_abs",
                         tau_levels = c(0, 0.5, 1), nlambda = 5,
                         nfolds = 3, parallel = FALSE),
               "tau_levels must lie strictly in \\(0, 1\\)")
})

test_that("cv.cpmnet rejects unsorted tau_levels", {
  dat <- sim_cpm()
  expect_error(cv.cpmnet(dat$x, dat$y, type = "pinball_abs",
                         tau_levels = c(0.5, 0.1, 0.9), nlambda = 5,
                         nfolds = 3, parallel = FALSE),
               "tau_levels must be sorted in increasing order")
})

test_that("cv.cpmnet rejects tau_weights not summing to 1", {
  dat <- sim_cpm()
  expect_error(cv.cpmnet(dat$x, dat$y, type = "pinball_abs",
                         tau_levels = c(0.25, 0.5, 0.75),
                         tau_weights = c(0.5, 0.5, 0.5),
                         nlambda = 5, nfolds = 3, parallel = FALSE),
               "tau_weights must sum to 1")
})

test_that("cv.cpmnet rejects brier_probs outside (0, 1)", {
  dat <- sim_cpm()
  expect_error(cv.cpmnet(dat$x, dat$y, type = "brier",
                         brier_probs = c(0.1, 0.5, 1.0),
                         nlambda = 5, nfolds = 3, parallel = FALSE),
               "brier_probs must lie strictly in \\(0, 1\\)")
})
