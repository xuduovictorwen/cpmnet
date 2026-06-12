#' Extract coefficients from a cpmnet fit
#'
#' Returns the slope coefficients (\code{beta}) and, optionally, the
#' ordinal intercepts (\code{alpha}) at one or all lambda values.
#'
#' @param object fitted \code{"cpmnet"} object.
#' @param s optional single lambda value at which to extract coefficients.
#'   If \code{NULL} (default), coefficients are returned for every lambda.
#'   If outside the fitted range, the nearest lambda is used with a warning.
#' @param include_alpha if \code{TRUE}, include the \code{k}-vector of
#'   ordinal intercepts. Default \code{FALSE} for glmnet-style slope-only
#'   output.
#' @param ... ignored; present for S3 compatibility.
#'
#' @return
#' When \code{s} is a scalar and \code{include_alpha = FALSE}: a named
#' numeric vector of length \code{p} (with predictor names when available).
#' When \code{s} is \code{NULL} and \code{include_alpha = FALSE}: a matrix
#' of dimension \code{p x length(object$lambda)}.
#' When \code{include_alpha = TRUE}: a list with components \code{alpha}
#' (intercepts) and \code{beta} (slopes) shaped as above.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{predict.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' fit <- cpmnet(x, y, family = "probit", nlambda = 10)
#' coef(fit, s = fit$lambda[5])
#' head(coef(fit))
#' }
#'
#' @export
coef.cpmnet <- function(object, s = NULL, include_alpha = FALSE, ...) {

  if (!inherits(object, "cpmnet")) {
    stop("object must inherit from class 'cpmnet'", call. = FALSE)
  }
  if (!is.logical(include_alpha) || length(include_alpha) != 1L ||
      is.na(include_alpha)) {
    stop("include_alpha must be a single TRUE/FALSE", call. = FALSE)
  }

  if (is.null(s)) {
    beta_out <- object$beta
    alpha_out <- object$alpha
    colnames(beta_out) <- paste0("s", seq_along(object$lambda))
    colnames(alpha_out) <- paste0("s", seq_along(object$lambda))
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
    idx <- which.min(abs(object$lambda - s))
    beta_out <- object$beta[, idx]
    alpha_out <- object$alpha[, idx]
    names(beta_out) <- rownames(object$beta)
    names(alpha_out) <- rownames(object$alpha)
  }

  if (include_alpha) {
    list(alpha = alpha_out, beta = beta_out)
  } else {
    beta_out
  }
}

#' Extract coefficients from a cv.cpmnet fit
#'
#' Returns slope coefficients at the lambda selected by a CV criterion.
#'
#' @param object fitted \code{"cv.cpmnet"} object.
#' @param s either a character string (\code{"lambda.min"} or
#'   \code{"lambda.1se"}) indicating which lambda to extract, or a
#'   numeric lambda value, or \code{NULL} to return the full path.
#'   For \code{type} in \code{"mean"}/\code{"median"} the names are
#'   \code{"lambda.min.mae"}, \code{"lambda.1se.mae"},
#'   \code{"lambda.min.mse"}, etc.; passing \code{"lambda.min"} maps to
#'   \code{"lambda.min.mae"} by default.
#' @param include_alpha if \code{TRUE}, also return the ordinal
#'   intercepts. Default \code{FALSE}.
#' @param ... ignored.
#'
#' @return As for \code{\link{coef.cpmnet}}, extracted at the selected
#'   lambda from \code{object$cpmnet.fit}.
#'
#' @seealso \code{\link{cv.cpmnet}}, \code{\link{coef.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5, type = "brier",
#'                    nlambda = 10, parallel = FALSE)
#' coef(cvfit, s = "lambda.min")
#' }
#'
#' @export
coef.cv.cpmnet <- function(object, s = "lambda.min",
                           include_alpha = FALSE, ...) {
  if (!inherits(object, "cv.cpmnet")) {
    stop("object must inherit from class 'cv.cpmnet'", call. = FALSE)
  }
  lambda_val <- .cv_lookup_lambda(object, s)
  coef.cpmnet(object$cpmnet.fit, s = lambda_val,
              include_alpha = include_alpha)
}

# Resolve a user-facing lambda spec (character or numeric, or NULL) into
# a numeric lambda value (or NULL for full-path). Handles the per-type
# name scheme used by cv.cpmnet.
#' @keywords internal
#' @noRd
.cv_lookup_lambda <- function(object, s) {
  if (is.null(s)) return(NULL)

  if (is.character(s)) {
    if (length(s) != 1L) {
      stop("s must be a single string or number", call. = FALSE)
    }
    # for mean/median, lambda.min / lambda.1se map to the MAE variant
    if (object$type %in% c("mean", "median")) {
      candidates <- list(
        "lambda.min"      = object$lambda.min.mae,
        "lambda.1se"      = object$lambda.1se.mae,
        "lambda.min.mae"  = object$lambda.min.mae,
        "lambda.1se.mae"  = object$lambda.1se.mae,
        "lambda.min.mse"  = object$lambda.min.mse,
        "lambda.1se.mse"  = object$lambda.1se.mse,
        "lambda.max.pR2"  = object$lambda.max.pR2,
        "lambda.1se.pR2"  = object$lambda.1se.pR2
      )
    } else {
      candidates <- list(
        "lambda.min" = object$lambda.min,
        "lambda.1se" = object$lambda.1se
      )
    }
    if (!(s %in% names(candidates))) {
      stop("Unknown lambda spec '", s, "'; valid options for type '",
           object$type, "': ",
           paste(names(candidates), collapse = ", "), call. = FALSE)
    }
    return(candidates[[s]])
  }

  if (is.numeric(s) && length(s) == 1L && is.finite(s)) return(s)

  stop("s must be NULL, a string, or a single number; got ",
       deparse(s), call. = FALSE)
}
