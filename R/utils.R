# replaces a trailing lambda_min of 0
.CPMNET_LAMBDA_MIN_FLOOR_RATIO <- 1e-8

# alpha_en 0 becomes 0.001 for lambda_max
.CPMNET_LAMBDA_MAX_ALPHA_FLOOR <- 0.001

#' @keywords internal
#' @noRd
validate_cpmnet_inputs <- function(x, y, family, alpha_en, weights, lambda,
                                   standardize = TRUE, adaptive = FALSE) {

  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("x must be a numeric matrix; got storage.mode '", storage.mode(x), "'",
         call. = FALSE)
  }
  n_bad_x <- sum(!is.finite(x))
  if (n_bad_x > 0) {
    stop("x must contain no NA/NaN/Inf; found ", n_bad_x, " non-finite value(s)",
         call. = FALSE)
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop("x must be non-empty; got ", nrow(x), " x ", ncol(x), " matrix",
         call. = FALSE)
  }

  if (!is.numeric(y)) {
    stop("y must be a numeric vector; got class '", paste(class(y), collapse = "/"),
         "'", call. = FALSE)
  }
  if (length(y) == 0L) stop("y must be non-empty", call. = FALSE)
  n_bad_y <- sum(!is.finite(y))
  if (n_bad_y > 0) {
    stop("y must contain no NA/NaN/Inf; found ", n_bad_y, " non-finite value(s)",
         call. = FALSE)
  }

  if (nrow(x) != length(y)) {
    stop("length(y) = ", length(y), " must equal nrow(x) = ", nrow(x),
         call. = FALSE)
  }

  if (!is.numeric(alpha_en) || length(alpha_en) != 1L || !is.finite(alpha_en)) {
    stop("alpha_en must be a single finite numeric value; got ",
         deparse(alpha_en), call. = FALSE)
  }
  if (alpha_en < 0 || alpha_en > 1) {
    stop("alpha_en must be in [0, 1]; got ", alpha_en, call. = FALSE)
  }

  valid_families <- names(rms::probabilityFamilies)
  if (!is.character(family) || length(family) != 1L) {
    stop("family must be a single character string; got ",
         deparse(family), call. = FALSE)
  }
  if (!family %in% valid_families) {
    stop("family must be one of: ", paste(valid_families, collapse = ", "),
         "; got '", family, "'", call. = FALSE)
  }

  if (!is.null(weights)) {
    if (!is.numeric(weights)) {
      stop("weights must be numeric; got class '",
           paste(class(weights), collapse = "/"), "'", call. = FALSE)
    }
    if (length(weights) != length(y)) {
      stop("length(weights) = ", length(weights),
           " must equal length(y) = ", length(y), call. = FALSE)
    }
    if (any(!is.finite(weights))) {
      stop("weights must contain no NA/NaN/Inf", call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative; found ", sum(weights < 0),
           " negative value(s)", call. = FALSE)
    }
    if (all(weights == 0)) {
      stop("weights must not all be zero", call. = FALSE)
    }
  }

  if (!is.null(lambda)) {
    if (!is.numeric(lambda)) {
      stop("lambda must be numeric; got class '",
           paste(class(lambda), collapse = "/"), "'", call. = FALSE)
    }
    if (any(!is.finite(lambda))) {
      stop("lambda must contain no NA/NaN/Inf", call. = FALSE)
    }
    if (any(lambda < 0)) {
      stop("lambda must be non-negative; found ", sum(lambda < 0),
           " negative value(s)", call. = FALSE)
    }
    if (length(lambda) >= 2L && is.unsorted(rev(lambda), strictly = FALSE)) {
      warning("lambda is not in decreasing order; sorting", call. = FALSE)
      lambda <- sort(lambda, decreasing = TRUE)
    }
  }

  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be a single TRUE/FALSE; got ",
         deparse(standardize), call. = FALSE)
  }
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("adaptive must be a single TRUE/FALSE; got ",
         deparse(adaptive), call. = FALSE)
  }

  list(x = x, y = y, weights = weights, lambda = lambda)
}

# drop zero weight rows, code y as categories 0..k, and set the intercepts of the
# model with no predictors, in the rms orientation
#' @keywords internal
#' @noRd
cpm_prep_y <- function(x, y, family, weights) {
  wt <- if (is.null(weights)) rep(1.0, length(y)) else weights
  if (any(wt == 0)) {
    active <- wt > 0
    x <- x[active, , drop = FALSE]
    y <- y[active]
    wt <- wt[active]
  }
  # rms::recode2integer is not exported
  recoded <- rms:::recode2integer(y, precision = 7)
  y_int <- recoded$y - 1L
  k <- max(y_int)
  if (k < 1L) {
    stop("y takes a single value (to 7 significant digits); a cumulative ",
         "probability model needs at least two", call. = FALSE)
  }
  sumwty <- tapply(wt, y_int, sum)
  cum_probs <- rev(cumsum(rev(sumwty)))[2:(k + 1)] / sum(sumwty)
  finverse <- eval(rms::probabilityFamilies[[family]][2])
  list(x = x, y = y_int, wt = wt, k = k, y_mapping = recoded$ylevels,
       link = match(family, names(rms::probabilityFamilies)),
       alpha_init = finverse(cum_probs))
}

#' @keywords internal
#' @noRd
cpm_data_prep <- function(x, y, family, weights, standardize,
                          x_means = NULL, x_sds = NULL) {

  prep <- cpm_prep_y(x, y, family, weights)
  x <- prep$x
  wt <- prep$wt
  y_int <- prep$y
  k <- prep$k
  n <- length(y_int)
  p <- ncol(x)

  if (standardize) {
    if (!is.null(x_means) && !is.null(x_sds)) {
      x <- t((t(x) - x_means) / x_sds)
    } else {
      is_weighted <- !all(wt == wt[1])
      if (is_weighted) {
        sum_wt <- sum(wt)
        x_means <- colSums(x * wt) / sum_wt
        x <- t(t(x) - x_means)
        x_vars <- colSums(wt * x * x) / sum_wt
      } else {
        x_means <- colMeans(x)
        x <- t(t(x) - x_means)
        x_vars <- colMeans(x * x)
      }
      x_sds <- sqrt(x_vars)
      x <- t(t(x) / x_sds)
    }
  } else {
    if (is.null(x_means)) x_means <- rep(0, p)
    if (is.null(x_sds)) x_sds <- rep(1, p)
  }

  link <- prep$link
  storage.mode(x) <- "double"
  storage.mode(y_int) <- "integer"
  storage.mode(wt) <- "double"
  storage.mode(k) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(link) <- "integer"

  list(x = x, y = y_int, wt = wt, n = n, p = p, k = k, link = link,
       alpha_init = prep$alpha_init, beta_init = rep(0, p),
       x_means = x_means, x_sds = x_sds, y_mapping = prep$y_mapping)
}

# grid floor, the fit separates below it
#' @keywords internal
#' @noRd
default_lambda_min_ratio <- function(n, p, family) {
  if (family == "cauchit") {
    if (n < 1.5 * p) 0.1 else if (n < 4 * p) 0.02 else 1e-4
  } else {
    if (n < 1.5 * p) 0.02 else 1e-4
  }
}

#' @keywords internal
#' @noRd
create_lambda_path <- function(lambda_max, lambda_min, nlambda) {
  if (lambda_min == 0) {
    return(c(exp(seq(log(lambda_max),
                     log(lambda_max * .CPMNET_LAMBDA_MIN_FLOOR_RATIO),
                     length.out = nlambda - 1)), 0))
  }
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}

# the smallest lambda at which every penalized coefficient is zero
#' @keywords internal
#' @noRd
compute_lambda_max <- function(data_prep, alpha_en, penalty_factor) {
  res <- .Fortran("cpm_lambda_max_grad",
                  data_prep$n, data_prep$k, data_prep$p,
                  data_prep$x, data_prep$y, data_prep$wt, data_prep$link,
                  data_prep$alpha_init,
                  u = numeric(data_prep$p), salloc = integer(1),
                  PACKAGE = "cpmnet")
  if (res$salloc != 0) stop("Fortran allocation failed in compute_lambda_max",
                            call. = FALSE)
  pen <- penalty_factor > 0
  max(abs(res$u[pen]) / penalty_factor[pen]) /
    max(alpha_en, .CPMNET_LAMBDA_MAX_ALPHA_FLOOR)
}

# after a Fortran path: stop on an allocation failure, keep the lambdas fitted before
# a singular intercept Hessian or a non-finite estimate, and report non-convergence
#' @keywords internal
#' @noRd
.path_status <- function(res, lambda, max_iter, verbose, routine) {
  if (res$salloc != 0 && res$salloc < 1000) {
    stop("Fortran allocation failed in ", routine, " (salloc = ", res$salloc, ")",
         call. = FALSE)
  }
  nl <- length(lambda)
  ok <- if (res$salloc >= 1000) which(res$niter < 0)[1] - 1L else nl
  keep <- seq_len(ok)
  bad <- apply(res$alpha_mat[, keep, drop = FALSE], 2, function(v) any(!is.finite(v))) |
         apply(res$beta_mat[, keep, drop = FALSE], 2, function(v) any(!is.finite(v)))
  if (any(bad)) ok <- which(bad)[1] - 1L
  if (ok == 0L) {
    if (res$salloc >= 1000) {
      stop("The intercept Hessian is singular at the first lambda (intercept ",
           res$salloc - 1000L, "). This usually means a nearly separable fit or ",
           "near-collinear intercepts. Try a larger lambda_min_ratio or remove ",
           "near-constant predictors.", call. = FALSE)
    }
    stop("Numerical issues at all lambda values. Try increasing lambda_min_ratio.",
         call. = FALSE)
  }
  if (ok < nl) {
    if (verbose) {
      why <- if (res$salloc >= 1000 && !any(bad)) "The intercept Hessian became singular"
             else "Numerical issues arose"
      warning(sprintf("%s at lambda %d of %d. Path truncated to %d lambdas.",
                      why, ok + 1L, nl, ok), call. = FALSE)
    }
    keep <- seq_len(ok)
    res$alpha_mat <- res$alpha_mat[, keep, drop = FALSE]
    res$beta_mat  <- res$beta_mat[, keep, drop = FALSE]
    res$logL      <- res$logL[keep]
    res$converged <- res$converged[keep]
    res$niter     <- res$niter[keep]
    lambda <- lambda[keep]
  }
  # a capped fit returns coefficients with no error, so the flag is the only signal
  if (verbose && any(!res$converged)) {
    warning(sprintf(
      "%d of %d lambda values hit max_iter = %d without converging; their coefficients are returned but not converged. See the `converged` field.",
      sum(!res$converged), length(res$converged), as.integer(max_iter)),
      call. = FALSE)
  }
  list(res = res, lambda = lambda)
}

# the coordinate visiting order is drawn inside the Fortran from a seed taken from
# R's stream, so set.seed pins a whole fit; the cpmnet.shuffle_seed option pins one
.cpmnet_shuffle_seed <- function() {
  s <- getOption("cpmnet.shuffle_seed", NULL)
  if (!is.null(s)) return(as.integer(s))
  as.integer(sample.int(.Machine$integer.max, 1L))
}
