# internal constants
# used in create_lambda_path when lambda_min == 0, we replace the trailing
# zero with a very small positive value for log-spaced stability
.CPMNET_LAMBDA_MIN_FLOOR_RATIO <- 1e-8

# used in compute_lambda_max to avoid division by zero when alpha_en is 0
# (pure ridge); matches the glmnet convention of treating alpha = 0 as
# alpha = 0.001 for lambda_max computation only
.CPMNET_LAMBDA_MAX_ALPHA_FLOOR <- 0.001

#' @keywords internal
#' @noRd
validate_cpmnet_inputs <- function(x, y, family, alpha_en, weights, lambda,
                                   standardize = TRUE, adaptive = FALSE) {

  # coerce x to matrix if data.frame/vector, then check numeric/finite
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

  # y numeric and finite
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

  # dimensional compatibility
  if (nrow(x) != length(y)) {
    stop("length(y) = ", length(y), " must equal nrow(x) = ", nrow(x),
         call. = FALSE)
  }

  # alpha_en in [0, 1]
  if (!is.numeric(alpha_en) || length(alpha_en) != 1L || !is.finite(alpha_en)) {
    stop("alpha_en must be a single finite numeric value; got ",
         deparse(alpha_en), call. = FALSE)
  }
  if (alpha_en < 0 || alpha_en > 1) {
    stop("alpha_en must be in [0, 1]; got ", alpha_en, call. = FALSE)
  }

  # family
  valid_families <- names(rms::probabilityFamilies)
  if (!is.character(family) || length(family) != 1L) {
    stop("family must be a single character string; got ",
         deparse(family), call. = FALSE)
  }
  if (!family %in% valid_families) {
    stop("family must be one of: ", paste(valid_families, collapse = ", "),
         "; got '", family, "'", call. = FALSE)
  }

  # weights
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

  # lambda
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

  # standardize and adaptive single logicals
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

#' @keywords internal
#' @noRd
cpm_data_prep <- function(x, y, family, weights, standardize,
                          x_means = NULL, x_sds = NULL) {

  wt <- if (is.null(weights)) rep(1.0, length(y)) else weights
  if (any(wt == 0)) {
    active <- wt > 0
    x <- x[active, , drop = FALSE]
    y <- y[active]
    wt <- wt[active]
  }

  # NOTE: rms::recode2integer is not exported. We depend on it for accurate
  # mapping of tied y values to integer ranks, and track rms releases to
  # detect breakage. See utils.R :: recode2integer in rms source.
  recode_result <- rms:::recode2integer(y, precision = 7)
  y_mapping <- recode_result$ylevels
  y_int <- recode_result$y - 1L

  n <- length(y_int)
  p <- ncol(x)
  k <- max(y_int)

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

  link <- match(family, names(rms::probabilityFamilies))
  sumwty <- tapply(wt, y_int, sum)
  cum_probs <- rev(cumsum(rev(sumwty)))[2:(k + 1)] / sum(sumwty)
  finverse <- eval(rms::probabilityFamilies[[family]][2])
  alpha_init <- finverse(cum_probs)
  if (any(!is.finite(alpha_init))) alpha_init <- seq(-2, 2, length.out = k)

  storage.mode(x) <- "double"
  storage.mode(y_int) <- "integer"
  storage.mode(wt) <- "double"
  storage.mode(k) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(link) <- "integer"

  list(x = x, y = y_int, wt = wt, n = n, p = p, k = k, link = link,
       alpha_init = alpha_init, beta_init = rep(0, p),
       x_means = x_means, x_sds = x_sds, y_mapping = y_mapping)
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

#' @keywords internal
#' @noRd
compute_lambda_max <- function(data_prep, alpha_en = 1) {
  fortran_fit <- .Fortran("ormll_beta_ccd_enet",
                          data_prep$n, data_prep$k, data_prep$p,
                          data_prep$x, data_prep$y, data_prep$wt,
                          link = data_prep$link,
                          data_prep$alpha_init, rep(0.0, data_prep$p),
                          lambda = 0.0, alpha_en = alpha_en,
                          logL = numeric(1),
                          u = numeric(data_prep$p),
                          hb = numeric(data_prep$p),
                          salloc = integer(1),
                          PACKAGE = "cpmnet")
  if (fortran_fit$salloc != 0) stop("Fortran allocation failed in compute_lambda_max",
                                    call. = FALSE)
  alpha_for_max <- max(alpha_en, .CPMNET_LAMBDA_MAX_ALPHA_FLOOR)
  max(abs(fortran_fit$u)) / alpha_for_max
}
