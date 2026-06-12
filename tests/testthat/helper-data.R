# Shared test fixtures.
# Using a seeded simulator keeps tests deterministic across sessions while
# leaving the global RNG state untouched for the surrounding test file.

sim_cpm <- function(n = 80, p = 5, sd = 1.0, seed = 42,
                    beta_true = NULL) {
  local_seed <- seed
  state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  set.seed(local_seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(beta_true)) {
    beta_true <- c(1, -1, rep(0, max(p - 2, 0)))[seq_len(p)]
  }
  y <- as.numeric(x %*% beta_true + stats::rnorm(n, sd = sd))
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv)
  }
  list(x = x, y = y, beta_true = beta_true, n = n, p = p)
}

# Small fit used by many tests: fast to compute and small enough that
# .Fortran calls stay cheap.
.small_fit <- function(family = "probit", nlambda = 10, alpha_en = 0.5,
                       n = 80, p = 5) {
  dat <- sim_cpm(n = n, p = p)
  cpmnet(dat$x, dat$y, family = family, nlambda = nlambda,
         alpha_en = alpha_en, verbose = FALSE)
}
