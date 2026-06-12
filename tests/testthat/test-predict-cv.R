# predict.cv.cpmnet: lambda string lookup, numeric s, NULL s, per-metric
# variants, error paths, delegation to predict.cpmnet.

.cv_pred_dat <- sim_cpm(n = 80, p = 4, seed = 300)

.cv_pred_fit <- function(type = "brier") {
  cv.cpmnet(.cv_pred_dat$x, .cv_pred_dat$y, type = type,
            nfolds = 3, nlambda = 8, parallel = FALSE,
            family = "probit")
}

test_that("predict.cv.cpmnet with s = 'lambda.min' returns 1-col matrix", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(cv, newx = newx, s = "lambda.min")
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, 1L))
})

test_that("predict.cv.cpmnet with s = 'lambda.1se' returns 1-col matrix", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(cv, newx = newx, s = "lambda.1se")
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, 1L))
})

test_that("predict.cv.cpmnet default s is 'lambda.min'", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  p_default <- predict(cv, newx = newx)
  p_min <- predict(cv, newx = newx, s = "lambda.min")
  expect_equal(p_default, p_min)
})

test_that("predict.cv.cpmnet with numeric s returns matrix at nearest lambda", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(cv, newx = newx, s = cv$lambda[3])
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, 1L))
})

test_that("predict.cv.cpmnet with s = NULL returns full path", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(cv, newx = newx, s = NULL)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, length(cv$lambda)))
})

test_that("predict.cv.cpmnet type = 'mean' differs from type = 'median' on skew", {
  set.seed(301)
  n <- 100; p <- 3
  x <- matrix(rnorm(n * p), n, p)
  y <- exp(as.numeric(x %*% c(1, -1, 0) + rnorm(n)))
  cv <- cv.cpmnet(x, y, family = "probit", type = "brier",
                  nfolds = 3, nlambda = 8, parallel = FALSE)
  newx <- matrix(rnorm(5 * p), 5, p)
  p_med  <- predict(cv, newx = newx, s = "lambda.min", type = "median")
  p_mean <- predict(cv, newx = newx, s = "lambda.min", type = "mean")
  expect_false(isTRUE(all.equal(p_med, p_mean)))
})

test_that("predict.cv.cpmnet resolves per-metric lambda names for type = 'mean'", {
  cv <- .cv_pred_fit("mean")
  newx <- matrix(rnorm(3 * 4), 3, 4)
  p_mae <- predict(cv, newx = newx, s = "lambda.min.mae")
  p_mse <- predict(cv, newx = newx, s = "lambda.min.mse")
  expect_true(is.matrix(p_mae))
  expect_true(is.matrix(p_mse))
  expect_equal(dim(p_mae), c(3L, 1L))
  expect_equal(dim(p_mse), c(3L, 1L))
})

test_that("predict.cv.cpmnet with type = 'mean' resolves lambda.max.pR2", {
  cv <- .cv_pred_fit("mean")
  newx <- matrix(rnorm(3 * 4), 3, 4)
  pred <- predict(cv, newx = newx, s = "lambda.max.pR2")
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(3L, 1L))
})

test_that("predict.cv.cpmnet errors on non-cv.cpmnet object", {
  expect_error(
    predict.cv.cpmnet(list(x = 1), newx = matrix(0, 1, 1)),
    "object must inherit from class 'cv.cpmnet'"
  )
})

test_that("predict.cv.cpmnet errors on unknown lambda string", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_error(
    predict(cv, newx = newx, s = "lambda.bogus"),
    "Unknown lambda spec"
  )
})

test_that("predict.cv.cpmnet matches manual predict.cpmnet delegation", {
  cv <- .cv_pred_fit("brier")
  newx <- matrix(rnorm(5 * 4), 5, 4)
  p_cv     <- predict(cv, newx = newx, s = "lambda.min")
  p_manual <- predict(cv$cpmnet.fit, newx = newx, s = cv$lambda.min,
                      type = "median")
  expect_equal(p_cv, p_manual)
})
