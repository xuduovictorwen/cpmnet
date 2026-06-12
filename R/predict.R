# Validate cpmnet object, coerce newx, and select lambda indices.
# Returns a list with newx (matrix), lambda_indices, and extracted
# alpha_mat / beta_mat / y_mapping / link_type ready for Fortran.
#' @keywords internal
#' @noRd
.prepare_predict_args <- function(object, newx, s = NULL) {

  if (!inherits(object, "cpmnet")) {
    stop("object must inherit from class 'cpmnet'; got class '",
         paste(class(object), collapse = "/"), "'", call. = FALSE)
  }

  # coerce newx
  if (is.vector(newx) && !is.matrix(newx)) newx <- matrix(newx, nrow = 1)
  if (!is.matrix(newx)) newx <- as.matrix(newx)
  if (!is.numeric(newx)) {
    stop("newx must be numeric; got storage.mode '", storage.mode(newx), "'",
         call. = FALSE)
  }
  if (any(!is.finite(newx))) {
    stop("newx must contain no NA/NaN/Inf", call. = FALSE)
  }
  expected_p <- nrow(object$beta)
  if (ncol(newx) != expected_p) {
    stop("ncol(newx) = ", ncol(newx),
         " must equal the number of predictors in the fit (",
         expected_p, ")", call. = FALSE)
  }

  # select lambda indices
  if (is.null(s)) {
    lambda_indices <- seq_along(object$lambda)
  } else {
    if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
      stop("s must be a single non-negative finite number; got ",
           deparse(s), call. = FALSE)
    }
    lam_rng <- range(object$lambda)
    if (s < lam_rng[1] || s > lam_rng[2]) {
      warning("s = ", s, " is outside the fitted lambda range [",
              lam_rng[1], ", ", lam_rng[2], "]; using nearest lambda",
              call. = FALSE)
    }
    lambda_indices <- which.min(abs(object$lambda - s))
  }

  alpha_mat <- object$alpha[, lambda_indices, drop = FALSE]
  beta_mat  <- object$beta[, lambda_indices, drop = FALSE]
  storage.mode(newx) <- "double"
  storage.mode(alpha_mat) <- "double"
  storage.mode(beta_mat) <- "double"
  y_mapping <- object$data_prep$y_mapping
  storage.mode(y_mapping) <- "double"

  list(
    newx = newx,
    n_obs = nrow(newx),
    p = ncol(newx),
    k = object$data_prep$k,
    n_lambda = length(lambda_indices),
    lambda_indices = lambda_indices,
    alpha_mat = alpha_mat,
    beta_mat = beta_mat,
    link_type = object$data_prep$link,
    y_mapping = y_mapping
  )
}

#' Predict from a cpmnet fit
#'
#' Compute conditional mean or median predictions from a fitted
#' \code{cpmnet} object at one or all lambda values.
#'
#' @param object fitted \code{"cpmnet"} object.
#' @param newx predictor matrix for new observations. A numeric vector is
#'   interpreted as a single observation (1 x p).
#' @param type either \code{"median"} (default) or \code{"mean"}.
#' @param s optional single lambda value at which to predict. If
#'   \code{NULL}, predictions are returned for every lambda in
#'   \code{object$lambda}. If outside the fitted range, the nearest
#'   lambda is used with a warning.
#' @param tol tolerance for the median interpolation at CDF = 0.5.
#' @param ... ignored; present for S3 compatibility.
#'
#' @return
#' A matrix of predictions with \code{nrow(newx)} rows and
#' \code{length(s)} (or \code{length(object$lambda)} when \code{s} is
#' \code{NULL}) columns. Column names are \code{"s<idx>=<lambda>"}.
#' Row names are copied from \code{newx} when present. Use
#' \code{drop()} or \code{as.vector()} to recover a plain vector.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100; p <- 5
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - x[, 2] + rnorm(n)
#' fit <- cpmnet(x, y, family = "probit", nlambda = 10)
#' # all lambdas, median prediction, first 3 observations
#' predict(fit, newx = x[1:3, ])
#' # single lambda = lambda.min placeholder
#' predict(fit, newx = x[1:3, ], type = "mean", s = fit$lambda[5])
#' }
#'
#' @export
predict.cpmnet <- function(object, newx, type = c("median", "mean"),
                           s = NULL, tol = 1e-10, ...) {
  type <- match.arg(type)
  pred_type <- if (type == "median") 1L else 2L
  args <- .prepare_predict_args(object, newx, s)
  result <- .Fortran("cpm_predict",
                     n = as.integer(args$n_obs), p = as.integer(args$p),
                     k = as.integer(args$k),
                     newx = args$newx, alpha = args$alpha_mat,
                     beta = args$beta_mat,
                     y_mapping = args$y_mapping,
                     link_type = as.integer(args$link_type),
                     pred_type = pred_type,
                     predictions = matrix(0.0, args$n_obs, args$n_lambda),
                     n_lambda = as.integer(args$n_lambda),
                     tol = as.double(tol))
  predictions <- result$predictions
  colnames(predictions) <- paste0("s", args$lambda_indices, "=",
                                  round(object$lambda[args$lambda_indices], 4))
  if (!is.null(rownames(args$newx))) rownames(predictions) <- rownames(args$newx)
  predictions
}

#' Predict from a cv.cpmnet fit
#'
#' Convenience wrapper that resolves \code{s} against the CV-selected
#' lambdas (\code{"lambda.min"}, \code{"lambda.1se"}, per-metric variants
#' for \code{type} in \code{"mean"}/\code{"median"}, a numeric value, or
#' \code{NULL} for the full path) and forwards to
#' \code{\link{predict.cpmnet}} on \code{object$cpmnet.fit}.
#'
#' @param object fitted \code{"cv.cpmnet"} object.
#' @param newx predictor matrix for new observations.
#' @param s either a character string indicating the CV-selected lambda
#'   (\code{"lambda.min"}, \code{"lambda.1se"}, or a per-metric variant
#'   like \code{"lambda.min.mae"}), a numeric lambda value, or
#'   \code{NULL} to predict at every lambda. Default \code{"lambda.min"}.
#' @param type either \code{"median"} (default) or \code{"mean"}.
#' @param tol tolerance passed to \code{\link{predict.cpmnet}}.
#' @param ... ignored; present for S3 compatibility.
#'
#' @return As for \code{\link{predict.cpmnet}}: an
#'   \code{nrow(newx) x n_selected_lambda} matrix.
#'
#' @seealso \code{\link{predict.cpmnet}}, \code{\link{coef.cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] - x[, 2] + rnorm(100)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5,
#'                    type = "brier", nlambda = 10, parallel = FALSE)
#' predict(cvfit, newx = x[1:5, ], s = "lambda.min")
#' predict(cvfit, newx = x[1:5, ], s = "lambda.1se", type = "mean")
#' }
#'
#' @export
predict.cv.cpmnet <- function(object, newx, s = "lambda.min",
                              type = c("median", "mean"),
                              tol = 1e-10, ...) {
  if (!inherits(object, "cv.cpmnet")) {
    stop("object must inherit from class 'cv.cpmnet'", call. = FALSE)
  }
  type <- match.arg(type)
  lambda_val <- .cv_lookup_lambda(object, s)
  predict.cpmnet(object$cpmnet.fit, newx = newx, type = type,
                 s = lambda_val, tol = tol)
}

# Quantile prediction via Fortran.
# Returns array of dimension n_obs x n_lambda x n_tau of predicted
# quantiles at the specified tau_levels.
#' @keywords internal
#' @noRd
predict_quantile_cpmnet <- function(object, newx, s = NULL, tau_levels,
                                    tol = 1e-10) {
  # tau_levels validation identical to cv.cpmnet
  if (!is.numeric(tau_levels) || length(tau_levels) == 0L) {
    stop("tau_levels must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(tau_levels))) {
    stop("tau_levels must contain no NA/NaN/Inf", call. = FALSE)
  }
  if (any(tau_levels <= 0) || any(tau_levels >= 1)) {
    stop("tau_levels must lie strictly in (0, 1)", call. = FALSE)
  }
  if (is.unsorted(tau_levels)) {
    stop("tau_levels must be sorted in increasing order", call. = FALSE)
  }

  args <- .prepare_predict_args(object, newx, s)
  storage.mode(tau_levels) <- "double"
  n_tau <- length(tau_levels)
  result <- .Fortran("cpm_predict_quantile",
                     n = as.integer(args$n_obs), p = as.integer(args$p),
                     k = as.integer(args$k),
                     newx = args$newx, alpha = args$alpha_mat,
                     beta = args$beta_mat,
                     y_mapping = args$y_mapping,
                     link_type = as.integer(args$link_type),
                     tau_levels = tau_levels,
                     n_tau = as.integer(n_tau),
                     predictions = array(0.0,
                                         dim = c(args$n_obs, args$n_lambda, n_tau)),
                     n_lambda = as.integer(args$n_lambda),
                     tol = as.double(tol))
  result$predictions
}

# Per-observation CDF for PSR.
# Returns n_obs x n_lambda matrix of CDF values F(c_i|X=x_i) evaluated
# at the supplied y_test. Splits by ordered support position (not
# numeric midpoint) to break ties at the boundary.
#' @keywords internal
#' @noRd
predict_psr_cpmnet <- function(object, newx, y_test, s = NULL) {

  args <- .prepare_predict_args(object, newx, s)

  if (!is.numeric(y_test) || length(y_test) != args$n_obs) {
    stop("y_test must be a numeric vector of length nrow(newx) = ",
         args$n_obs, "; got length ", length(y_test), call. = FALSE)
  }
  if (any(!is.finite(y_test))) {
    stop("y_test must contain no NA/NaN/Inf", call. = FALSE)
  }

  k <- args$k
  y_map <- object$data_prep$y_mapping
  # split by ordered support position, not numeric midpoint
  j <- findInterval(y_test, y_map)
  mid_j <- floor(k / 2L)
  obs_idx <- integer(args$n_obs)
  for (i in seq_len(args$n_obs)) {
    if (j[i] == 0L) {
      obs_idx[i] <- 1L
    } else if (j[i] >= k + 1L) {
      obs_idx[i] <- k + 1L
    } else if (y_test[i] == y_map[j[i]]) {
      obs_idx[i] <- j[i]
    } else if (j[i] <= mid_j) {
      obs_idx[i] <- min(j[i] + 1L, k + 1L)
    } else {
      obs_idx[i] <- j[i]
    }
  }
  result <- .Fortran("cpm_predict_psr",
                     n = as.integer(args$n_obs), p = as.integer(args$p),
                     k = as.integer(k),
                     newx = args$newx, alpha = args$alpha_mat,
                     beta = args$beta_mat,
                     link_type = as.integer(args$link_type),
                     obs_idx = as.integer(obs_idx),
                     cdf_out = matrix(0.0, args$n_obs, args$n_lambda),
                     n_lambda = as.integer(args$n_lambda))
  result$cdf_out
}

# CDF prediction via Fortran.
# Returns array of dimension n_obs x n_lambda x n_thresh of CDF values
# F(c|X=x_i) evaluated at the supplied thresholds.
#' @keywords internal
#' @noRd
predict_cdf_cpmnet <- function(object, newx, s = NULL, thresholds) {

  if (!is.numeric(thresholds) || length(thresholds) == 0L) {
    stop("thresholds must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(thresholds))) {
    stop("thresholds must contain no NA/NaN/Inf", call. = FALSE)
  }

  args <- .prepare_predict_args(object, newx, s)
  storage.mode(thresholds) <- "double"
  n_thresh <- length(thresholds)
  result <- .Fortran("cpm_predict_cdf",
                     n = as.integer(args$n_obs), p = as.integer(args$p),
                     k = as.integer(args$k),
                     newx = args$newx, alpha = args$alpha_mat,
                     beta = args$beta_mat,
                     y_mapping = args$y_mapping,
                     link_type = as.integer(args$link_type),
                     thresholds = thresholds,
                     n_thresh = as.integer(n_thresh),
                     cdf_out = array(0.0,
                                     dim = c(args$n_obs, args$n_lambda, n_thresh)),
                     n_lambda = as.integer(args$n_lambda))
  result$cdf_out
}
