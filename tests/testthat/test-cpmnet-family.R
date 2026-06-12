# One fit per supported link family. We check class + finite output rather
# than coefficient values, which vary by link.

valid_families <- c("logistic", "probit", "loglog", "cloglog", "cauchit")

for (fam in valid_families) {
  local({
    fam_local <- fam
    test_that(sprintf("family '%s' fits and returns finite alpha/beta",
                      fam_local), {
      dat <- sim_cpm(n = 80, p = 4)
      fit <- cpmnet(dat$x, dat$y, family = fam_local, nlambda = 10,
                    verbose = FALSE)
      expect_s3_class(fit, "cpmnet")
      expect_true(all(is.finite(fit$alpha)))
      expect_true(all(is.finite(fit$beta)))
    })
  })
}
