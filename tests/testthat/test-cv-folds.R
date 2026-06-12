# Fold handling: user-supplied foldid, various nfolds, seed determinism.

test_that("user-supplied foldid is used exactly", {
  dat <- sim_cpm(n = 30, p = 4, seed = 50)
  fid <- rep_len(1:3, 30)
  cv <- cv.cpmnet(dat$x, dat$y, foldid = fid, type = "brier",
                  nlambda = 6, nfolds = 3, parallel = FALSE,
                  family = "probit")
  expect_equal(cv$nfolds_successful, 3L)
})

test_that("nfolds = 2 is accepted", {
  dat <- sim_cpm(n = 40, p = 4, seed = 51)
  cv <- cv.cpmnet(dat$x, dat$y, nfolds = 2, nlambda = 5,
                  type = "brier", parallel = FALSE, family = "probit")
  expect_s3_class(cv, "cv.cpmnet")
  expect_equal(cv$nfolds_successful, 2L)
})

test_that("nfolds = 5 is accepted", {
  dat <- sim_cpm(n = 60, p = 4, seed = 52)
  cv <- cv.cpmnet(dat$x, dat$y, nfolds = 5, nlambda = 5,
                  type = "brier", parallel = FALSE, family = "probit")
  expect_equal(cv$nfolds_successful, 5L)
})

test_that("fixed seed gives deterministic CV result (sequential)", {
  dat <- sim_cpm(n = 40, p = 4, seed = 60)
  set.seed(999)
  cv1 <- cv.cpmnet(dat$x, dat$y, nfolds = 3, nlambda = 6,
                   type = "brier", parallel = FALSE, family = "probit")
  set.seed(999)
  cv2 <- cv.cpmnet(dat$x, dat$y, nfolds = 3, nlambda = 6,
                   type = "brier", parallel = FALSE, family = "probit")
  expect_equal(cv1$cvm, cv2$cvm)
  expect_equal(cv1$lambda.min, cv2$lambda.min)
})

test_that("foldid with non-sequential max works", {
  dat <- sim_cpm(n = 30, p = 4, seed = 53)
  fid <- rep_len(1:3, 30)
  cv <- cv.cpmnet(dat$x, dat$y, foldid = fid, nlambda = 5,
                  type = "brier", parallel = FALSE, nfolds = 3,
                  family = "probit")
  expect_equal(max(fid), 3L)
  expect_equal(cv$nfolds_successful, 3L)
})

test_that("duplicate-seed CV with user foldid is fully deterministic", {
  dat <- sim_cpm(n = 40, p = 4, seed = 54)
  fid <- rep_len(1:4, 40)
  cv1 <- cv.cpmnet(dat$x, dat$y, foldid = fid, nlambda = 6,
                   type = "brier", parallel = FALSE, nfolds = 4,
                   family = "probit")
  cv2 <- cv.cpmnet(dat$x, dat$y, foldid = fid, nlambda = 6,
                   type = "brier", parallel = FALSE, nfolds = 4,
                   family = "probit")
  expect_equal(cv1$cvm, cv2$cvm)
})

test_that("foldid as double coerces to integer if exact", {
  dat <- sim_cpm(n = 30, p = 4, seed = 55)
  fid <- as.double(rep_len(1:3, 30))
  cv <- cv.cpmnet(dat$x, dat$y, foldid = fid, nlambda = 5,
                  type = "brier", parallel = FALSE, nfolds = 3,
                  family = "probit")
  expect_s3_class(cv, "cv.cpmnet")
})

test_that("cv.cpmnet forwards dots arguments to cpmnet", {
  dat <- sim_cpm(n = 40, p = 4, seed = 56)
  cv <- cv.cpmnet(dat$x, dat$y, nfolds = 3, nlambda = 5,
                  type = "brier", parallel = FALSE,
                  family = "probit", alpha_en = 0.1)
  expect_equal(cv$cpmnet.fit$alpha_en, 0.1)
})
