#' Fit a penalized (regularized) cumulative probability model
#'
#' Fits a cumulative probability model (CPM) with elastic net penalty
#' along a penalization (regularization) path using coordinate descent. The core loop is
#' implemented in Fortran; R drives setup, standardization, and
#' post-processing.
#'
#' @param x predictor matrix with \code{n} rows and \code{p} columns.
#'   Must be numeric with no NA/NaN/Inf. Coerced to matrix if a data frame.
#' @param y response vector of length \code{n}. Numeric, treated as
#'   ordinal; ties are handled by ranking.
#' @param family link function name. One of \code{"logistic"},
#'   \code{"probit"}, \code{"loglog"}, \code{"cloglog"}, \code{"cauchit"}
#'   (or any name in \code{names(rms::probabilityFamilies)}).
#' @param lambda optional user-supplied lambda sequence (decreasing).
#'   If \code{NULL}, a path is computed automatically from \code{lambda_max}
#'   down to \code{lambda_max * lambda_min_ratio}.
#' @param alpha_en elastic net mixing parameter in \code{[0, 1]}.
#'   \code{1} = lasso, \code{0} = ridge. See Note about naming.
#' @param nlambda number of lambda values when \code{lambda} is \code{NULL}.
#'   Must be an integer >= 2.
#' @param lambda_min_ratio ratio of \code{lambda.min} to \code{lambda.max}.
#'   Defaults by link and by \code{n} against \code{p}, since the fit separates
#'   below the floor. For \code{"cauchit"}, \code{0.1} if \code{n < 1.5p},
#'   \code{0.02} if \code{n < 4p}, \code{1e-4} otherwise. For the other four
#'   links, \code{0.02} if \code{n < 1.5p}, \code{1e-4} otherwise.
#'   Must be positive.
#' @param lambda_max optional override for the maximum lambda. If
#'   \code{NULL}, computed from the gradient at the intercept-only fit.
#' @param weights observation weights (non-negative) of length \code{n}.
#'   Defaults to uniform weights. Observations with weight 0 are dropped.
#' @param penalty_factor per-predictor penalty multipliers of length \code{p}.
#'   Non-negative; not all zero. Defaults to \code{rep(1, p)}.
#' @param adaptive if \code{TRUE} and \code{penalty_factor} is \code{NULL},
#'   compute adaptive lasso weights from an unpenalized \code{rms::orm.fit}.
#' @param standardize if \code{TRUE}, standardize predictors to unit
#'   variance before fitting and unstandardize coefficients at return.
#' @param max_iter maximum outer iterations per lambda value.
#' @param tol convergence tolerance on max absolute parameter change.
#' @param nulldev null deviance. Computed from intercept-only
#'   \code{rms::orm.fit} if \code{NULL}.
#' @param x_means optional precomputed column means for standardization.
#' @param x_sds optional precomputed column standard deviations.
#' @param verbose if \code{TRUE} (default), emit a warning when the
#'   lambda path is truncated due to numerical issues.
#' @param ... additional arguments. The name \code{alpha} is accepted as
#'   a synonym for \code{alpha_en} for glmnet-style ergonomics; other
#'   names produce an error.
#'
#' @return An object of class \code{"cpmnet"}: a list with components
#' \describe{
#'   \item{alpha}{\code{k x nlambda} matrix of intercepts on the original
#'     x-scale (after un-standardization).}
#'   \item{beta}{\code{p x nlambda} matrix of slopes on the original
#'     x-scale.}
#'   \item{alpha_scaled, beta_scaled}{Same as \code{alpha}, \code{beta}
#'     but on the standardized scale used by the Fortran core.}
#'   \item{df}{Number of nonzero coefficients per lambda.}
#'   \item{dim}{\code{c(p, nlambda)}.}
#'   \item{lambda}{The lambda sequence used.}
#'   \item{dev.ratio}{\code{1 - logL / nulldev}, clipped at zero.}
#'   \item{nulldev, logL}{Null and model deviances.}
#'   \item{converged, niter}{Per-lambda convergence flags and iteration counts.}
#'   \item{nobs}{Number of observations after dropping zero-weight rows.}
#'   \item{data_prep}{Internal prep bundle (used by \code{predict.cpmnet}).}
#'   \item{alpha_en, standardize, call}{Echo of arguments and the call.}
#' }
#'
#' @note
#' The argument is named \code{alpha_en} rather than \code{alpha} because
#' \code{alpha} is used throughout the CPM literature and in this package
#' to denote the \code{k}-vector of ordinal intercepts stored on the
#' returned fit object. To ease transition from glmnet, \code{alpha} can
#' still be passed via \code{...} and will be forwarded to \code{alpha_en}.
#'
#' The returned object carries a copy of \code{data_prep}, which includes
#' the (standardized) design matrix. For memory-sensitive workflows,
#' set \code{object$data_prep$x <- NULL} after fitting if predictions are
#' no longer needed.
#'
#' @references
#' Liu, Q., Shepherd, B. E., Li, C., and Harrell, F. E. (2017).
#' Modeling continuous response variables using ordinal regression.
#' \emph{Statistics in Medicine}, 36(27), 4316-4335.
#'
#' Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization
#' paths for generalized linear models via coordinate descent.
#' \emph{Journal of Statistical Software}, 33(1), 1-22.
#'
#' @seealso \code{\link{cv.cpmnet}}, \code{\link{predict.cpmnet}},
#'   \code{\link{coef.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(1, -1, 0.5, rep(0, p - 3))
#' y <- x %*% beta_true + rnorm(n)
#' fit <- cpmnet(x, y, family = "probit", alpha_en = 0.5, nlambda = 20)
#' print(fit)
#' plot(fit, xvar = "lambda")
#' }
#'
#' @export
cpmnet <- function(x, y, family = "probit", lambda = NULL, alpha_en = 0.5,
                   nlambda = 100, lambda_min_ratio = NULL, lambda_max = NULL,
                   weights = NULL, penalty_factor = NULL, adaptive = FALSE,
                   standardize = TRUE, max_iter = 1000,
                   tol = 1e-7, nulldev = NULL,
                   x_means = NULL, x_sds = NULL, verbose = TRUE, ...) {

  # alpha is a synonym for alpha_en
  dots <- list(...)
  if ("alpha" %in% names(dots)) {
    if (!missing(alpha_en)) {
      stop("Both `alpha` and `alpha_en` supplied; use only one.",
           call. = FALSE)
    }
    alpha_en <- dots$alpha
    dots$alpha <- NULL
  }
  if (length(dots) > 0L) {
    stop("Unused arguments: ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }

  validated <- validate_cpmnet_inputs(x, y, family, alpha_en, weights, lambda,
                                      standardize, adaptive)
  lambda <- validated$lambda  # may have been sorted

  # scalar hyperparameter checks
  if (!is.numeric(nlambda) || length(nlambda) != 1L ||
      !is.finite(nlambda) || nlambda != as.integer(nlambda) || nlambda < 2L) {
    stop("nlambda must be a single integer >= 2; got ", deparse(nlambda),
         call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter != as.integer(max_iter) || max_iter < 1L) {
    stop("max_iter must be a single positive integer; got ", deparse(max_iter),
         call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a single positive number; got ", deparse(tol),
         call. = FALSE)
  }
  if (!is.null(lambda_min_ratio)) {
    if (!is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1L ||
        !is.finite(lambda_min_ratio) || lambda_min_ratio < 0) {
      stop("lambda_min_ratio must be a single non-negative number; got ",
           deparse(lambda_min_ratio), call. = FALSE)
    }
  }
  if (!is.null(lambda_max)) {
    if (!is.numeric(lambda_max) || length(lambda_max) != 1L ||
        !is.finite(lambda_max) || lambda_max <= 0) {
      stop("lambda_max must be a single positive finite number; got ",
           deparse(lambda_max), call. = FALSE)
    }
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be a single TRUE/FALSE; got ", deparse(verbose),
         call. = FALSE)
  }

  data_prep <- cpm_data_prep(validated$x, validated$y, family,
                             validated$weights, standardize, x_means, x_sds)

  # penalty_factor validation
  if (!is.null(penalty_factor)) {
    if (!is.numeric(penalty_factor)) {
      stop("penalty_factor must be numeric; got class '",
           paste(class(penalty_factor), collapse = "/"), "'", call. = FALSE)
    }
    if (length(penalty_factor) != data_prep$p) {
      stop("length(penalty_factor) = ", length(penalty_factor),
           " must equal p = ", data_prep$p, call. = FALSE)
    }
    if (any(!is.finite(penalty_factor))) {
      stop("penalty_factor must contain no NA/NaN/Inf", call. = FALSE)
    }
    if (any(penalty_factor < 0)) {
      stop("penalty_factor must be non-negative; found ",
           sum(penalty_factor < 0), " negative value(s)", call. = FALSE)
    }
    if (all(penalty_factor == 0)) {
      stop("penalty_factor must not be all zero", call. = FALSE)
    }
  }

  if (adaptive && is.null(penalty_factor)) {
    orm_fit <- rms::orm.fit(
      x = data_prep$x,
      y = data_prep$y_mapping[data_prep$y + 1],
      family = family, weights = data_prep$wt)
    orm_beta <- orm_fit$coefficients[-(1:data_prep$k)]
    penalty_factor <- 1.0 / abs(orm_beta)
  }

  if (is.null(penalty_factor)) penalty_factor <- rep(1.0, data_prep$p)
  storage.mode(penalty_factor) <- "double"

  if (is.null(lambda)) {
    if (is.null(lambda_max))
      lambda_max <- compute_lambda_max(data_prep, alpha_en)
    if (is.null(lambda_min_ratio))
      lambda_min_ratio <- default_lambda_min_ratio(data_prep$n, data_prep$p, family)
    if (lambda_min_ratio == 0 && data_prep$n <= data_prep$p)
      stop("lambda_min_ratio = 0 not allowed when n <= p", call. = FALSE)
    lambda_path <- create_lambda_path(
      lambda_max, lambda_max * lambda_min_ratio, nlambda)
  } else {
    lambda_path <- sort(lambda, decreasing = TRUE)
  }
  nlambda <- length(lambda_path)

  fit_result <- .Fortran("ormll_path_enet",
                         n = data_prep$n,
                         k = data_prep$k,
                         p = data_prep$p,
                         x = data_prep$x,
                         y = data_prep$y,
                         wt = data_prep$wt,
                         link = data_prep$link,
                         nlambda = as.integer(nlambda),
                         lambda_path = lambda_path,
                         alpha_en = alpha_en,
                         penalty_factor = penalty_factor,
                         alpha_init = data_prep$alpha_init,
                         beta_init = data_prep$beta_init,
                         max_iter_outer = as.integer(max_iter),
                         tol = tol,
                         alpha_mat = matrix(0.0, data_prep$k, nlambda),
                         beta_mat = matrix(0.0, data_prep$p, nlambda),
                         logL = numeric(nlambda),
                         converged = logical(nlambda),
                         niter = integer(nlambda),
                         salloc = integer(1),
                         PACKAGE = "cpmnet")

  if (fit_result$salloc != 0) {
    if (fit_result$salloc >= 1000) {
      stop("Singular Hessian in Fortran tridiagonal solve at alpha index ",
           fit_result$salloc - 1000L,
           ". This usually means degenerate y (near-collinear intercepts) or",
           " an ill-conditioned design. Try increasing lambda_min_ratio or",
           " removing near-constant predictors.", call. = FALSE)
    }
    stop("Fortran allocation failed in ormll_path_enet (salloc = ",
         fit_result$salloc, ")", call. = FALSE)
  }

  # truncate path at first NaN
  nan_lambda <- apply(fit_result$alpha_mat, 2, function(col) any(!is.finite(col))) |
                apply(fit_result$beta_mat, 2, function(col) any(!is.finite(col)))

  if (any(nan_lambda)) {
    first_nan <- which(nan_lambda)[1]
    valid_idx <- seq_len(first_nan - 1)
    if (length(valid_idx) == 0)
      stop("Numerical issues at all lambda values. Try increasing lambda_min_ratio.",
           call. = FALSE)
    if (verbose) {
      warning(sprintf("Numerical issues at %d lambda value(s). Path truncated to %d lambdas.",
                      sum(nan_lambda), length(valid_idx)), call. = FALSE)
    }
    fit_result$alpha_mat <- fit_result$alpha_mat[, valid_idx, drop = FALSE]
    fit_result$beta_mat  <- fit_result$beta_mat[, valid_idx, drop = FALSE]
    fit_result$logL      <- fit_result$logL[valid_idx]
    fit_result$converged <- fit_result$converged[valid_idx]
    fit_result$niter     <- fit_result$niter[valid_idx]
    lambda_path <- lambda_path[valid_idx]
    nlambda <- length(valid_idx)
  }

  df <- colSums(fit_result$beta_mat != 0)

  # unstandardize
  if (standardize) {
    beta_mat <- fit_result$beta_mat / data_prep$x_sds
    adjustment <- as.numeric(crossprod(data_prep$x_means, beta_mat))
    alpha_mat <- fit_result$alpha_mat -
      matrix(adjustment, nrow = data_prep$k, ncol = nlambda, byrow = TRUE)
  } else {
    beta_mat  <- fit_result$beta_mat
    alpha_mat <- fit_result$alpha_mat
  }

  if (is.null(nulldev)) {
    orm_null <- rms::orm.fit(
      x = NULL, y = data_prep$y_mapping[data_prep$y + 1],
      family = family, weights = data_prep$wt)
    nulldev <- orm_null$deviance
  }

  structure(
    list(
      alpha       = alpha_mat,
      beta        = beta_mat,
      alpha_scaled = fit_result$alpha_mat,
      beta_scaled  = fit_result$beta_mat,
      df          = df,
      dim         = c(data_prep$p, nlambda),
      lambda      = lambda_path,
      dev.ratio   = pmax(0, 1 - fit_result$logL / nulldev),
      nulldev     = nulldev,
      logL        = fit_result$logL,
      converged   = fit_result$converged,
      niter       = fit_result$niter,
      nobs        = data_prep$n,
      data_prep   = data_prep,
      alpha_en    = alpha_en,
      standardize = standardize,
      call        = match.call()
    ),
    class = "cpmnet"
  )
}
