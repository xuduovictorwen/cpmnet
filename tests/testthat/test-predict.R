# predict.cpmnet across return shapes and types.

.pred_fit <- function() .small_fit(nlambda = 10, n = 80, p = 4)

test_that("predict with s = NULL returns matrix nrow(newx) x nlambda", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(fit, newx = newx)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, length(fit$lambda)))
})

test_that("predict always returns matrix, even with single-row newx", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(4), 1, 4)
  pred <- predict(fit, newx = newx)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(1L, length(fit$lambda)))
})

test_that("predict with scalar s returns n x 1 matrix", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(fit, newx = newx, s = fit$lambda[3])
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(5L, 1L))
})

test_that("predict with scalar s and single-row newx returns 1 x 1 matrix", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(4), 1, 4)
  pred <- predict(fit, newx = newx, s = fit$lambda[3])
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(1L, 1L))
})

test_that("predict type = 'mean' returns numeric output", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  pred <- predict(fit, newx = newx, type = "mean")
  expect_true(all(is.finite(pred)))
})

test_that("predict type = 'median' returns numeric output", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  pred <- predict(fit, newx = newx, type = "median")
  expect_true(all(is.finite(pred)))
})

test_that("predict median and mean differ on non-symmetric data", {
  set.seed(200)
  n <- 100; p <- 3
  x <- matrix(rnorm(n * p), n, p)
  y <- exp(as.numeric(x %*% c(1, -1, 0) + rnorm(n)))  # skewed
  fit <- cpmnet(x, y, family = "probit", nlambda = 10, verbose = FALSE)
  newx <- matrix(rnorm(5 * p), 5, p)
  m1 <- predict(fit, newx = newx, type = "median", s = fit$lambda[5])
  m2 <- predict(fit, newx = newx, type = "mean", s = fit$lambda[5])
  expect_false(isTRUE(all.equal(m1, m2)))
})

test_that("newx as numeric vector treated as 1-row matrix", {
  fit <- .pred_fit()
  newx_vec <- as.numeric(rnorm(4))
  pred <- predict(fit, newx = newx_vec)
  expect_equal(dim(pred), c(1L, length(fit$lambda)))
})

test_that("predict errors on ncol(newx) != p", {
  fit <- .pred_fit()
  newx_bad <- matrix(rnorm(5 * 2), 5, 2)
  expect_error(predict(fit, newx = newx_bad),
               "ncol\\(newx\\) = 2 must equal the number of predictors")
})

test_that("predict errors on non-finite newx", {
  fit <- .pred_fit()
  newx_bad <- matrix(rnorm(5 * 4), 5, 4); newx_bad[1, 1] <- NA
  expect_error(predict(fit, newx = newx_bad),
               "newx must contain no NA/NaN/Inf")
})

test_that("predict warns when s is outside fitted lambda range", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_warning(
    predict(fit, newx = newx, s = 1e6),
    "outside the fitted lambda range"
  )
})

test_that("predict errors on negative s", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(3 * 4), 3, 4)
  expect_error(predict(fit, newx = newx, s = -0.1),
               "s must be a single non-negative finite number")
})

test_that("predict on non-cpmnet object errors", {
  expect_error(predict.cpmnet(list(x = 1), newx = matrix(0, 1, 1)),
               "object must inherit from class 'cpmnet'")
})

test_that("column names reflect the lambda index and value", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(5 * 4), 5, 4)
  pred <- predict(fit, newx = newx)
  expect_true(all(grepl("^s[0-9]+=", colnames(pred))))
})

test_that("predict preserves row names of newx", {
  fit <- .pred_fit()
  newx <- matrix(rnorm(5 * 4), 5, 4)
  rownames(newx) <- paste0("obs", 1:5)
  pred <- predict(fit, newx = newx)
  expect_equal(rownames(pred), rownames(newx))
})
