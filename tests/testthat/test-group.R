# grouped penalization. the point of each test is that a factor must behave as
# one object: selected or dropped whole, and unaffected by which level was
# dropped as the reference.

make_factor_data <- function(n = 240, seed = 3) {
  set.seed(seed)
  z <- sample(rep(1:4, each = n / 4))
  mu <- c(0, 0.8, 1.6, 2.4)
  xc <- matrix(rnorm(n * 4), n, 4)
  y <- exp(mu[z] + 0.5 * xc[, 1] + rnorm(n))
  list(x = cbind(model.matrix(~ factor(z))[, -1, drop = FALSE], xc),
       y = y, group = c(1, 1, 1, 2, 3, 4, 5), z = z, xc = xc)
}

test_that("group = NULL and an all-singleton group take the same route", {
  d <- make_factor_data()
  # the visiting order is drawn fresh per fit, so pin it to compare the two fits
  old <- options(cpmnet.shuffle_seed = 1L); on.exit(options(old))
  a <- cpmnet(d$x, d$y, family = "probit", alpha_en = 0.5, verbose = FALSE)
  b <- cpmnet(d$x, d$y, family = "probit", alpha_en = 0.5,
              group = seq_len(ncol(d$x)), verbose = FALSE)
  # identical, not merely close: an all-singleton grouping must not change the core
  expect_identical(a$beta, b$beta)
  expect_identical(a$lambda, b$lambda)
})

test_that("a grouped factor enters and leaves as one block", {
  d <- make_factor_data()
  f <- cpmnet(d$x, d$y, family = "probit", alpha_en = 1, group = d$group,
              verbose = FALSE)
  blk <- f$beta[1:3, , drop = FALSE]
  nz <- colSums(blk != 0)
  # every lambda has all three dummies in or all three out, never a subset
  expect_true(all(nz %in% c(0L, 3L)))
})

test_that("the factor's penalty does not depend on the reference level", {
  d <- make_factor_data()
  lam <- cpmnet(d$x, d$y, family = "probit", alpha_en = 1, group = d$group,
                verbose = FALSE)$lambda
  norms <- vapply(1:4, function(b) {
    levs <- c(b, setdiff(1:4, b))
    xb <- cbind(model.matrix(~ factor(d$z, levels = levs))[, -1, drop = FALSE], d$xc)
    f <- cpmnet(xb, d$y, family = "probit", alpha_en = 1, group = d$group,
                lambda = lam, max_iter = 200000L, tol = 1e-10, verbose = FALSE)
    il <- which(colSums(f$beta_tilde[1:3, , drop = FALSE] != 0) > 0)[1]
    sqrt(sum(f$beta_tilde[1:3, il]^2))
  }, numeric(1))
  expect_gt(norms[1], 0)
  # the whole claim of the method: four codings, one penalty
  expect_lt(max(norms) - min(norms), 1e-6)
})

test_that("group 0 leaves a column unpenalized", {
  d <- make_factor_data()
  g <- d$group; g[4] <- 0                     # first continuous column unpenalized
  f <- cpmnet(d$x, d$y, family = "probit", alpha_en = 1, group = g, verbose = FALSE)
  # an unpenalized column is never shrunk to zero, including at the largest lambda
  expect_true(all(f$beta[4, ] != 0))
})

test_that("the structural predict types report groups, not columns", {
  d <- make_factor_data()
  f <- cpmnet(d$x, d$y, family = "probit", alpha_en = 1, group = d$group,
              verbose = FALSE)
  ng <- predict(f, type = "ngroups")
  nv <- predict(f, type = "nvars")
  nm <- predict(f, type = "norm")
  expect_length(ng, length(f$lambda))
  expect_equal(nrow(nm), 5L)                  # five groups, not seven columns
  # a live factor counts once as a group and three times as columns
  expect_true(all(nv >= ng))
  expect_identical(unname(ng), unname(f$ngroups))
})

test_that("a grouped fit still works with the existing methods", {
  d <- make_factor_data()
  f <- cpmnet(d$x, d$y, family = "probit", alpha_en = 1, group = d$group,
              verbose = FALSE)
  expect_silent(cf <- coef(f, s = f$lambda[10]))
  pr <- predict(f, newx = d$x, s = f$lambda[10])
  expect_equal(nrow(pr), nrow(d$x))
  expect_output(print(f))
})

test_that("a bad group argument is rejected", {
  d <- make_factor_data()
  expect_error(cpmnet(d$x, d$y, family = "probit", group = c(1, 1), verbose = FALSE),
               "must equal ncol")
})
