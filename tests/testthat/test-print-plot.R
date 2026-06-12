# print / plot / summary / coef surface tests. These are mostly "does it
# run and return the right shape" checks since the rendered output is not
# part of the public contract.

test_that("print.cpmnet returns the object invisibly", {
  fit <- .small_fit()
  expect_invisible(invisible(capture.output(out <- print(fit))))
  # and in practice:
  out <- capture.output(ret <- print(fit))
  expect_identical(ret, fit)
})

test_that("print.cv.cpmnet runs for every type", {
  dat <- sim_cpm(n = 60, p = 4, seed = 600)
  for (ty in c("brier", "mean", "median", "pinball_abs", "psr_abs")) {
    cv <- cv.cpmnet(dat$x, dat$y, type = ty, nlambda = 5,
                    nfolds = 3, parallel = FALSE, family = "probit")
    expect_silent(capture.output(print(cv)))
  }
})

test_that("plot.cpmnet runs for xvar in {lambda, norm, dev}", {
  fit <- .small_fit(nlambda = 15)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (xv in c("lambda", "norm", "dev")) {
    expect_silent(plot(fit, xvar = xv))
  }
})

test_that("plot.cv.cpmnet runs for a Brier CV object", {
  dat <- sim_cpm(n = 60, p = 4, seed = 601)
  cv <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 6,
                  nfolds = 3, parallel = FALSE, family = "probit")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(cv))
})

test_that("plot.cv.cpmnet rejects bad metric with a helpful message", {
  dat <- sim_cpm(n = 60, p = 4, seed = 602)
  cv <- cv.cpmnet(dat$x, dat$y, type = "mean", nlambda = 5,
                  nfolds = 3, parallel = FALSE, family = "probit")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(cv, metric = "foo"),
               "metric for type 'mean' must be one of: mae, mse, pR2")
})

test_that("summary.cpmnet returns a data frame invisibly", {
  fit <- .small_fit()
  out <- capture.output(tbl <- summary(fit))
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("index", "lambda", "df", "dev.ratio") %in% names(tbl)))
})

test_that("summary.cv.cpmnet runs and returns a list", {
  dat <- sim_cpm(n = 60, p = 4, seed = 603)
  cv <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 5,
                  nfolds = 3, parallel = FALSE, family = "probit")
  out <- capture.output(lst <- summary(cv))
  expect_type(lst, "list")
  expect_true("lambda.min" %in% names(lst))
})

test_that("coef.cpmnet with s = NULL returns p x nlambda matrix", {
  fit <- .small_fit()
  co <- coef(fit)
  expect_equal(dim(co), c(nrow(fit$beta), length(fit$lambda)))
})

test_that("coef.cpmnet with scalar s returns named vector of length p", {
  fit <- .small_fit()
  co <- coef(fit, s = fit$lambda[3])
  expect_equal(length(co), nrow(fit$beta))
})

test_that("coef.cpmnet include_alpha = TRUE returns list(alpha, beta)", {
  fit <- .small_fit()
  out <- coef(fit, s = fit$lambda[3], include_alpha = TRUE)
  expect_type(out, "list")
  expect_named(out, c("alpha", "beta"))
  expect_equal(length(out$alpha), fit$data_prep$k)
  expect_equal(length(out$beta), nrow(fit$beta))
})

test_that("coef.cv.cpmnet with 'lambda.min' returns vector", {
  dat <- sim_cpm(n = 60, p = 4, seed = 604)
  cv <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 6,
                  nfolds = 3, parallel = FALSE, family = "probit")
  co <- coef(cv, s = "lambda.min")
  expect_equal(length(co), 4L)
})

test_that("coef.cv.cpmnet rejects unknown lambda spec", {
  dat <- sim_cpm(n = 60, p = 4, seed = 605)
  cv <- cv.cpmnet(dat$x, dat$y, type = "brier", nlambda = 5,
                  nfolds = 3, parallel = FALSE, family = "probit")
  expect_error(coef(cv, s = "lambda.foo"),
               "Unknown lambda spec")
})
