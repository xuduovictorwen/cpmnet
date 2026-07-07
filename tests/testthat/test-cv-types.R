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

test_that("loglik returns expected structure", {
  cv <- .run_cv("loglik")
  expect_s3_class(cv, "cv.cpmnet")
  expect_equal(cv$type, "loglik")
  expect_true(is.finite(cv$lambda.min))
  expect_true(is.finite(cv$lambda.1se))
  # cvm is a negative mean log of category masses in (0, 1], so >= 0
  expect_true(all(cv$cvm >= 0, na.rm = TRUE))
})

test_that("loglik category mass matches CDF differences and boundary clamps", {
  fit <- cpmnet(.cv_dat$x, .cv_dat$y, family = "probit", nlambda = 5)
  y_map <- fit$data_prep$y_mapping
  k <- fit$data_prep$k
  s <- fit$lambda[3]
  newx <- .cv_dat$x[1:2, , drop = FALSE]
  # y at interior training atoms: mass = F(u_j) - F(u_j-1)
  jsel <- c(3L, 5L)
  mass <- cpmnet:::predict_loglik_mass_cpmnet(fit, newx, y_map[jsel], s = s)
  cdf <- cpmnet:::predict_cdf_cpmnet(fit, newx, s = s,
                                     thresholds = y_map[1:k])
  manual <- c(cdf[1, 1, jsel[1]] - cdf[1, 1, jsel[1] - 1L],
              cdf[2, 1, jsel[2]] - cdf[2, 1, jsel[2] - 1L])
  expect_equal(as.numeric(mass), manual, tolerance = 1e-12)
  # below the support -> F(u_1); above -> 1 - F(u_k)
  y_out <- c(min(y_map) - 1, max(y_map) + 1)
  mass_out <- cpmnet:::predict_loglik_mass_cpmnet(fit, newx, y_out, s = s)
  expect_equal(mass_out[1, 1], cdf[1, 1, 1], tolerance = 1e-12)
  expect_equal(mass_out[2, 1], 1 - cdf[2, 1, k], tolerance = 1e-12)
  # masses are probabilities
  expect_true(all(mass > 0 & mass <= 1))
  expect_true(all(mass_out > 0 & mass_out <= 1))
})

test_that("cv.cpmnet accepts verbose passed through to cpmnet", {
  cv <- cv.cpmnet(.cv_dat$x, .cv_dat$y, nfolds = 3, nlambda = 8,
                  parallel = FALSE, family = "probit", verbose = FALSE)
  expect_s3_class(cv, "cv.cpmnet")
  expect_true(is.finite(cv$lambda.min))
})

test_that("cpmnet.fit is attached and is a valid cpmnet object", {
  cv <- .run_cv("brier")
  expect_s3_class(cv$cpmnet.fit, "cpmnet")
})

test_that("lambda.min is one of the fitted lambdas", {
  cv <- .run_cv("brier")
  expect_true(cv$lambda.min %in% cv$lambda)
})
