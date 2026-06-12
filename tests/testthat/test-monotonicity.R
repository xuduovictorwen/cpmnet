# CPM correctness: intercepts alpha_j must be strictly decreasing so that
# P(Y >= y_j | x) = F(alpha_j + x*beta) is non-increasing in j. The CCD
# update does not explicitly project onto this constraint; empirically it
# is preserved on clean data across all link families and sample sizes we
# tested (see plan item 60). This test guards against regressions.

test_that("fitted intercepts alpha are strictly decreasing for every family", {
  set.seed(99)
  families <- c("logistic", "probit", "loglog", "cloglog", "cauchit")
  for (fam in families) {
    for (n in c(60, 120)) {
      x <- matrix(rnorm(n * 4), n, 4)
      y <- as.numeric(x %*% c(1, -1, 0.5, 0) + rnorm(n))
      fit <- cpmnet(x, y, family = fam, nlambda = 10, verbose = FALSE)
      # every column of alpha (each lambda) must have all negative diffs
      worst_diff <- max(apply(fit$alpha, 2, function(a) max(diff(a))))
      expect_lt(worst_diff, 0, label = sprintf(
        "family = %s, n = %d: max(diff(alpha)) = %.4f (must be < 0)",
        fam, n, worst_diff
      ))
    }
  }
})

test_that("alpha_scaled (standardized scale) is also strictly decreasing", {
  set.seed(100)
  dat <- sim_cpm(n = 80, p = 4)
  fit <- cpmnet(dat$x, dat$y, family = "probit", nlambda = 10,
                verbose = FALSE)
  worst_diff <- max(apply(fit$alpha_scaled, 2, function(a) max(diff(a))))
  expect_lt(worst_diff, 0)
})
