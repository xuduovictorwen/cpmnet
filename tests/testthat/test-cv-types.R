# One cv.cpmnet run per type value. Shared dataset keeps timing tractable.

# Small shared dataset so each type runs quickly.
.cv_dat <- sim_cpm(n = 80, p = 4, seed = 100)

.cv_common <- list(
  nfolds = 3, nlambda = 8, parallel = FALSE, family = "probit"
)

.run_cv <- function(type) {
  do.call(cv.cpmnet,
          c(list(x = .cv_dat$x, y = .cv_dat$y, type = type),
            .cv_common))
}

test_that("brier (default) returns expected structure", {
  cv <- .run_cv("brier")
  expect_s3_class(cv, "cv.cpmnet")
  expect_equal(cv$type, "brier")
  expect_true(is.finite(cv$lambda.min))
  expect_true(is.finite(cv$lambda.1se))
  expect_true(all(c("cvm", "cvsd", "brier_probs", "brier_cuts") %in%
                  names(cv)))
})

test_that("brier default is selected without specifying type", {
  cv <- cv.cpmnet(.cv_dat$x, .cv_dat$y, nfolds = 3, nlambda = 8,
                  parallel = FALSE, family = "probit")
  expect_equal(cv$type, "brier")
})

test_that("mean type produces MAE/MSE/pR2 structure", {
  cv <- .run_cv("mean")
  expect_equal(cv$type, "mean")
  expect_true(all(c("cvm.mae", "cvm.mse", "cvm.pR2") %in% names(cv)))
  expect_true(is.finite(cv$lambda.min.mae))
  expect_true(is.finite(cv$lambda.min.mse))
})

test_that("median type produces MAE/MSE (pR2 all NA)", {
  cv <- .run_cv("median")
  expect_equal(cv$type, "median")
  expect_true(is.finite(cv$lambda.min.mae))
  # pR2 is not computed for median, so all NA
  expect_true(all(is.na(cv$cvm.pR2)))
})

test_that("pinball_abs returns pinball structure", {
  cv <- .run_cv("pinball_abs")
  expect_equal(cv$type, "pinball_abs")
  expect_true(is.finite(cv$lambda.min))
  expect_true("tau_levels" %in% names(cv))
})

test_that("pinball_sq returns pinball structure", {
  cv <- .run_cv("pinball_sq")
  expect_equal(cv$type, "pinball_sq")
  expect_true(is.finite(cv$lambda.min))
})

test_that("psr_abs returns psr structure", {
  cv <- .run_cv("psr_abs")
  expect_equal(cv$type, "psr_abs")
  expect_true(is.finite(cv$lambda.min))
  # PSR values lie in [0, 1]
  expect_true(all(cv$cvm >= 0 & cv$cvm <= 1, na.rm = TRUE))
})

test_that("psr_sq returns psr structure", {
  cv <- .run_cv("psr_sq")
  expect_equal(cv$type, "psr_sq")
  expect_true(is.finite(cv$lambda.min))
  expect_true(all(cv$cvm >= 0 & cv$cvm <= 1, na.rm = TRUE))
})

test_that("cpmnet.fit is attached and is a valid cpmnet object", {
  cv <- .run_cv("brier")
  expect_s3_class(cv$cpmnet.fit, "cpmnet")
})

test_that("lambda.min is one of the fitted lambdas", {
  cv <- .run_cv("brier")
  expect_true(cv$lambda.min %in% cv$lambda)
})
