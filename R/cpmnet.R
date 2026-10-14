#' Fit a penalized (regularized) cumulative probability model
#'
#' Fits a cumulative probability model (CPM) with elastic net penalty
#' along a penalization (regularization) path. The intercepts are updated
#' together by a damped Newton step on their tridiagonal Hessian, and the
#' coefficients by coordinate or group descent. The core loop is
#' implemented in Fortran; R drives setup, standardization, and
#' post-processing.
#'
#' @param x predictor matrix with \code{n} rows and \code{p} columns.
#'   Must be numeric with no NA/NaN/Inf. Coerced to matrix if a data frame.
#'   A column that takes a single value gets coefficient 0 and is left out
#'   of the fit.
#' @param y response vector of length \code{n}. Numeric, treated as
#'   ordinal; ties are handled by ranking. Must take at least two values.
#' @param family link function name. One of \code{"logistic"},
#'   \code{"probit"}, \code{"loglog"}, \code{"cloglog"}, \code{"cauchit"}
#'   (or any name in \code{names(rms::probabilityFamilies)}).
#' @param group optional length-\code{p} vector assigning each column of
#'   \code{x} to a group, so that a factor is penalized as one unit. Accepts a
#'   factor, character or integer vector. Columns labelled \code{0} are left
#'   unpenalized. \code{NULL}, the default, and any grouping in which every
#'   penalized group holds a single column both give an elementwise penalty. A
#'   grouping with any block of two or more columns centers and orthonormalizes
#'   each block, so the penalty is on the factor's fitted contribution and does
#'   not depend on which level was dropped.
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
#'   Both cores center and scale the predictors with the weights.
#' @param penalty_factor penalty multipliers, non-negative and not all zero:
#'   one per column of \code{x}, or one per penalized group when \code{group}
#'   holds a block of two or more columns. Defaults to 1 for each.
#' @param adaptive if \code{TRUE} and \code{penalty_factor} is \code{NULL},
#'   compute adaptive lasso weights from an unpenalized \code{rms::orm.fit}.
#' @param core which solver to use when the penalty is elementwise. Both
#'   minimize the same objective. \code{"group"}, the default, is faster on
#'   logistic, probit and cauchit and slower on the two extreme value links.
#'   \code{"elementwise"} updates one coefficient at a time. The two agree to
#'   about 1e-7 on the four log concave links; on cauchit they can differ more,
#'   because its objective is not convex and the two can reach different local
#'   minima. A grouping with a block of two or more columns always uses the
#'   grouped core.
#' @param standardize if \code{TRUE}, standardize predictors to unit
#'   variance before fitting and unstandardize coefficients at return.
#' @param max_iter maximum outer iterations per lambda. Default 100000,
#'   matching \pkg{glmnet}. A warning is issued if any lambda reaches it.
#' @param tol convergence tolerance on max absolute parameter change.
#' @param nulldev null deviance. Computed from intercept-only
#'   \code{rms::orm.fit} if \code{NULL}.
#' @param x_means optional precomputed column means for standardization.
#' @param x_sds optional precomputed column standard deviations.
#' @param verbose if \code{TRUE} (default), warn when the path is truncated,
#'   when a lambda does not converge, or when a constant column is dropped.
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
#'   \item{constant_columns}{Indices of the columns of \code{x} that were
#'     constant and left out of the fit, when there are any.}
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
cpmnet <- function(x, y, family = "probit", group = NULL,
                   lambda = NULL, alpha_en = 0.5,
                   nlambda = 100, lambda_min_ratio = NULL, lambda_max = NULL,
                   weights = NULL, penalty_factor = NULL, adaptive = FALSE,
                   standardize = TRUE, max_iter = 100000,
                   tol = 1e-7, nulldev = NULL,
                   x_means = NULL, x_sds = NULL, core = c("group", "elementwise"),
                   verbose = TRUE, ...) {

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

  core <- match.arg(core)

  # one factor per column without a grouping, one per penalized group with one
  if (!is.null(penalty_factor)) {
    if (!is.numeric(penalty_factor)) {
      stop("penalty_factor must be numeric; got class '",
           paste(class(penalty_factor), collapse = "/"), "'", call. = FALSE)
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

  # adaptive weights, no standardization and supplied means or sds exist only in
  # the elementwise core; a block of two or more columns needs the grouped one
  elementwise_only <- adaptive || !standardize ||
    !is.null(x_means) || !is.null(x_sds)
  p <- ncol(validated$x)
  gspec <- if (is.null(group)) NULL else normalize_group(group, p)
  real_group <- !is.null(gspec) &&
    (!gspec$all_singleton || gspec$n_unpenalized > 0L)

  if (!is.null(penalty_factor)) {
    np <- if (real_group) length(gspec$labels) else p
    if (length(penalty_factor) != np) {
      if (real_group) {
        stop("penalty_factor must have length = number of groups (", np, ")",
             call. = FALSE)
      }
      stop("length(penalty_factor) = ", length(penalty_factor),
           " must equal p = ", np, call. = FALSE)
    }
  }

  if (real_group && elementwise_only) {
    stop("A grouping cannot be combined with adaptive, standardize = FALSE, ",
         "x_means or x_sds; those belong to the elementwise core.",
         call. = FALSE)
  }

  # a constant column has nothing to fit; it gets coefficient 0 and stays out
  rows <- if (is.null(validated$weights)) TRUE else validated$weights > 0
  varies <- apply(validated$x[rows, , drop = FALSE], 2, function(v) any(v != v[1]))
  if (!any(varies)) {
    stop("every column of x is constant, so there is nothing to fit", call. = FALSE)
  }
  x_fit <- validated$x
  gspec_fit <- gspec
  pf_fit <- penalty_factor
  if (!all(varies)) {
    if (verbose) {
      warning(sprintf("%d constant column(s) of x get coefficient 0: %s",
                      sum(!varies), paste(utils::head(which(!varies), 10), collapse = ", ")),
              call. = FALSE)
    }
    x_fit <- validated$x[, varies, drop = FALSE]
    if (!is.null(x_means)) x_means <- x_means[varies]
    if (!is.null(x_sds)) x_sds <- x_sds[varies]
    if (!is.null(gspec)) gspec_fit <- normalize_group(group[varies], sum(varies))
    if (!is.null(penalty_factor)) {
      pf_fit <- if (real_group) penalty_factor[match(gspec_fit$labels, gspec$labels)]
                else penalty_factor[varies]
      if (all(pf_fit == 0)) {
        stop("every penalized column of x is constant", call. = FALSE)
      }
    }
  }

  if (real_group || (core == "group" && !elementwise_only)) {
    gid <- if (is.null(gspec_fit)) seq_len(ncol(x_fit)) else gspec_fit$gid
    attr(gid, "labels") <- if (is.null(gspec_fit)) NULL else gspec_fit$labels
    fit <- cpmnet_group_path(
      x = x_fit, y = validated$y, gid = gid, family = family,
      alpha_en = alpha_en, lambda = validated$lambda, nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio, lambda_max = lambda_max,
      penalty_factor = pf_fit, weights = validated$weights,
      max_iter = as.integer(max_iter), tol = tol, nulldev = nulldev,
      verbose = verbose, mc = match.call())
  } else {
    fit <- cpmnet_elementwise_path(
      x = x_fit, y = validated$y, family = family, alpha_en = alpha_en,
      lambda = validated$lambda, nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio, lambda_max = lambda_max,
      weights = validated$weights, penalty_factor = pf_fit, adaptive = adaptive,
      standardize = standardize, max_iter = as.integer(max_iter), tol = tol,
      nulldev = nulldev, x_means = x_means, x_sds = x_sds, verbose = verbose,
      mc = match.call())
  }

  if (!all(varies)) {
    nl <- length(fit$lambda)
    beta <- matrix(0, p, nl)
    beta[varies, ] <- fit$beta
    fit$beta <- beta
    if (is.null(fit$group)) {
      bs <- matrix(0, p, nl)
      bs[varies, ] <- fit$beta_scaled
      fit$beta_scaled <- bs
    } else {
      gid <- if (is.null(gspec)) seq_len(p) else gspec$gid
      attr(gid, "labels") <- if (is.null(gspec)) NULL else gspec$labels
      fit$group <- gid
      fit$group_labels <- attr(gid, "labels")
    }
    fit$dim <- c(p, nl)
    fit$constant_columns <- which(!varies)
  }
  fit
}

# the elementwise path, one coefficient at a time
#' @keywords internal
#' @noRd
cpmnet_elementwise_path <- function(x, y, family, alpha_en, lambda, nlambda,
                                    lambda_min_ratio, lambda_max, weights,
                                    penalty_factor, adaptive, standardize,
                                    max_iter, tol, nulldev, x_means, x_sds,
                                    verbose, mc) {
  data_prep <- cpm_data_prep(x, y, family, weights, standardize, x_means, x_sds)

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

  first_null <- FALSE
  if (is.null(lambda)) {
    # at a computed lambda_max every coefficient is zero; a floored alpha_en or an
    # unpenalized column breaks that
    first_null <- is.null(lambda_max) && alpha_en >= .CPMNET_LAMBDA_MAX_ALPHA_FLOOR &&
      all(penalty_factor > 0)
    if (is.null(lambda_max))
      lambda_max <- compute_lambda_max(data_prep, alpha_en, penalty_factor)
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
  # rounding at lambda_max can leave a coefficient near 1e-15, so the first fit takes a
  # far larger penalty and keeps every coefficient exactly zero, as glmnet does
  lambda_fit <- lambda_path
  if (first_null) lambda_fit[1] <- 1e10 * lambda_path[1]

  res <- .Fortran("ormll_path_enet",
                  n = data_prep$n,
                  k = data_prep$k,
                  p = data_prep$p,
                  x = data_prep$x,
                  y = data_prep$y,
                  wt = data_prep$wt,
                  link = data_prep$link,
                  nlambda = as.integer(nlambda),
                  lambda_path = lambda_fit,
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
                  shuffle_seed = .cpmnet_shuffle_seed(),
                  PACKAGE = "cpmnet")
  st <- .path_status(res, lambda_path, max_iter, verbose, "ormll_path_enet")
  res <- st$res
  lambda_path <- st$lambda
  nlambda <- length(lambda_path)

  df <- colSums(res$beta_mat != 0)

  if (standardize) {
    beta_mat <- res$beta_mat / data_prep$x_sds
    adjustment <- as.numeric(crossprod(data_prep$x_means, beta_mat))
    alpha_mat <- res$alpha_mat -
      matrix(adjustment, nrow = data_prep$k, ncol = nlambda, byrow = TRUE)
  } else {
    beta_mat  <- res$beta_mat
    alpha_mat <- res$alpha_mat
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
      alpha_scaled = res$alpha_mat,
      beta_scaled  = res$beta_mat,
      df          = df,
      dim         = c(data_prep$p, nlambda),
      lambda      = lambda_path,
      dev.ratio   = pmax(0, 1 - res$logL / nulldev),
      nulldev     = nulldev,
      logL        = res$logL,
      converged   = res$converged,
      niter       = res$niter,
      nobs        = data_prep$n,
      data_prep   = data_prep,
      alpha_en    = alpha_en,
      standardize = standardize,
      call        = mc
    ),
    class = "cpmnet"
  )
}
