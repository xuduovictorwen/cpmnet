#' Cross-validation for cpmnet
#'
#' Performs k-fold cross-validation for a \code{cpmnet} model and selects
#' lambda according to one of five conceptual criteria: mean, median,
#' pinball, Brier, or held-out log-likelihood. The default, \code{"brier"},
#' averages the squared Brier score across a user-specified grid of
#' response percentiles and is the criterion recommended in the
#' accompanying methods paper.
#'
#' @param x predictor matrix. See \code{\link{cpmnet}} for requirements.
#' @param y response vector. See \code{\link{cpmnet}}.
#' @param nfolds number of folds. Must satisfy \code{2 <= nfolds <= n}.
#' @param foldid optional integer vector of length \code{n} assigning each
#'   observation to a fold. When supplied, values must form
#'   \code{1:max(foldid)} with no gaps.
#' @param type cross-validation loss. One of \code{"brier"} (default),
#'   \code{"mean"}, \code{"median"}, \code{"pinball_abs"},
#'   \code{"pinball_sq"}, \code{"psr_abs"}, \code{"psr_sq"},
#'   \code{"loglik"}.
#'   These map to five conceptual criteria emphasized in the paper:
#'   mean and median regression loss (with MAE/MSE/pR2 reported),
#'   pinball quantile loss (absolute or squared variant), the
#'   Brier score, and the held-out log-likelihood. The \code{"psr_*"}
#'   options remain in the package for completeness but are not part
#'   of the paper's primary set. \code{"loglik"} selects the lambda
#'   that maximizes the mean held-out log-likelihood of the observed
#'   outcome categories (the criterion \code{ordinalNet} calls
#'   \code{cvLoglik}); \code{cvm} stores its negative so that
#'   \code{lambda.min} keeps its usual minimizing meaning. A held-out
#'   outcome outside the training support is assigned the boundary
#'   category, following the same convention as the \code{"psr_*"}
#'   criteria; category masses are floored at \code{1e-12} before the
#'   log.
#' @param parallel if \code{TRUE} (default), use \code{parallel::mclapply}.
#'   Ignored on Windows with a warning.
#' @param ncores number of cores. Defaults to
#'   \code{min(nfolds, parallel::detectCores() - 1)}.
#' @param tau_levels quantile levels used by \code{type = "pinball_*"}.
#'   Must be a sorted numeric vector in \code{(0, 1)}. Default is
#'   \code{seq(0.1, 0.9, by = 0.1)}.
#' @param tau_weights weights for each tau level; non-negative, sum to 1.
#'   Default is uniform.
#' @param brier_probs percentile levels defining Brier thresholds when
#'   \code{type = "brier"}. Must be sorted, in \code{(0, 1)}. Default is
#'   \code{seq(0.1, 0.9, by = 0.1)}.
#' @param brier_weights weights for each Brier threshold; non-negative,
#'   sum to 1. Default is uniform.
#' @param ... arguments passed to \code{\link{cpmnet}}.
#'
#' @return An object of class \code{"cv.cpmnet"}. Structure depends on
#'   \code{type}; all variants carry \code{lambda}, \code{lambda.min} /
#'   \code{lambda.1se} (or \code{.mae}, \code{.mse}, \code{.pR2} suffixes
#'   for \code{type} in \code{"mean"}/\code{"median"}), \code{cpmnet.fit}
#'   (full fit on all data), and \code{index}.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{predict.cpmnet}},
#'   \code{\link{coef.cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 120; p <- 5
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - x[, 2] + rnorm(n)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5, type = "brier",
#'                    nlambda = 20, parallel = FALSE)
#' cvfit$lambda.min
#' }
#'
#' @export
cv.cpmnet <- function(x, y, nfolds = 10, foldid = NULL,
                      type = c("brier", "mean", "median", "pinball_abs",
                               "pinball_sq", "psr_abs", "psr_sq", "loglik"),
                      parallel = TRUE, ncores = NULL,
                      tau_levels = seq(0.1, 0.9, by = 0.1),
                      tau_weights = NULL,
                      brier_probs = seq(0.1, 0.9, by = 0.1),
                      brier_weights = NULL, ...) {

  type <- match.arg(type)

  # dimensional baseline
  if (is.vector(x) && !is.matrix(x)) x <- as.matrix(x)
  n <- length(y)

  # nfolds
  if (!is.numeric(nfolds) || length(nfolds) != 1L ||
      !is.finite(nfolds) || nfolds != as.integer(nfolds)) {
    stop("nfolds must be a single integer; got ", deparse(nfolds),
         call. = FALSE)
  }
  if (nfolds < 2L) {
    stop("nfolds must be at least 2; got ", nfolds, call. = FALSE)
  }
  if (nfolds > n) {
    stop("nfolds = ", nfolds, " exceeds n = ", n, call. = FALSE)
  }

  # foldid
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length.out = n))
  } else {
    if (!is.numeric(foldid) || length(foldid) != n) {
      stop("foldid must be an integer vector of length n = ", n,
           "; got length ", length(foldid), call. = FALSE)
    }
    if (any(!is.finite(foldid)) || any(foldid != as.integer(foldid))) {
      stop("foldid must contain only integer values", call. = FALSE)
    }
    foldid <- as.integer(foldid)
    uniq <- sort(unique(foldid))
    if (min(uniq) < 1L || !identical(uniq, seq_len(max(uniq)))) {
      stop("foldid values must form the set 1:max(foldid) with no gaps; got ",
           paste(uniq, collapse = ", "), call. = FALSE)
    }
  }

  # ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) || length(ncores) != 1L ||
        !is.finite(ncores) || ncores != as.integer(ncores) || ncores < 1L) {
      stop("ncores must be a single positive integer; got ", deparse(ncores),
           call. = FALSE)
    }
  }

  # parallel
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be a single TRUE/FALSE; got ", deparse(parallel),
         call. = FALSE)
  }

  # probability grids and weights
  tau_weights <- .check_prob_grid(tau_levels, tau_weights,
                                  "tau_levels", "tau_weights")
  n_tau <- length(tau_levels)
  brier_weights <- .check_prob_grid(brier_probs, brier_weights,
                                    "brier_probs", "brier_weights")
  n_brier <- length(brier_probs)

  dots <- list(...)
  fit_full <- do.call(cpmnet, c(list(x = x, y = y), dots))
  lambda_path <- fit_full$lambda
  nlambda <- length(lambda_path)
  # brier thresholds before splitting
  brier_cuts <- as.numeric(quantile(y, probs = brier_probs))

  process_fold <- function(fold) {
    val_idx <- which(foldid == fold)
    train_idx <- which(foldid != fold)
    n_val <- length(val_idx)
    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_val <- x[val_idx, , drop = FALSE]
    y_val <- y[val_idx]
    # fold fits are silent
    fold_args <- c(list(x = x_train, y = y_train, lambda = lambda_path), dots)
    if (is.null(fold_args$verbose)) fold_args$verbose <- FALSE
    fit_fold <- do.call(cpmnet, fold_args)
    nlambda_fold <- length(fit_fold$lambda)
    if (type %in% c("mean", "median")) {
      # returns n_val x nlambda_fold
      pred <- predict.cpmnet(fit_fold, newx = x_val, type = type, s = NULL)
      mae <- rep(NA_real_, nlambda)
      mse <- rep(NA_real_, nlambda)
      pR2 <- rep(NA_real_, nlambda)
      mae[1:nlambda_fold] <- colMeans(abs(pred - y_val))
      mse[1:nlambda_fold] <- colMeans((pred - y_val)^2)
      if (type == "mean") {
        for (j in 1:nlambda_fold) {
          yhat <- pred[, j]
          PAm <- PAmeasures::pam.nlm(y_val, yhat)
          pR2[j] <- as.numeric(PAm[1]) * as.numeric(PAm[2])
        }
      }
      list(mae = mae, mse = mse, pR2 = pR2)
    } else if (type %in% c("pinball_abs", "pinball_sq")) {
      qpred <- predict_quantile_cpmnet(fit_fold, newx = x_val,
                                        s = NULL, tau_levels = tau_levels)
      pinball <- rep(NA_real_, nlambda)
      for (j in 1:nlambda_fold) {
        pb_tau <- numeric(n_tau)
        for (t in 1:n_tau) {
          tau <- tau_levels[t]
          qhat <- qpred[, j, t]
          resid <- y_val - qhat
          if (type == "pinball_abs") {
            pb_tau[t] <- mean(ifelse(resid >= 0, tau * abs(resid), (1 - tau) * abs(resid)))
          } else {
            pb_tau[t] <- mean(ifelse(resid >= 0, tau * resid^2, (1 - tau) * resid^2))
          }
        }
        pinball[j] <- sum(tau_weights * pb_tau)
      }
      list(pinball = pinball)
    } else if (type == "brier") {
      cdf_pred <- predict_cdf_cpmnet(fit_fold, newx = x_val,
                                      s = NULL, thresholds = brier_cuts)
      brier <- rep(NA_real_, nlambda)
      for (j in 1:nlambda_fold) {
        bs_c <- numeric(n_brier)
        for (t in 1:n_brier) {
          fhat <- cdf_pred[, j, t]
          indicator <- as.numeric(y_val <= brier_cuts[t])
          bs_c[t] <- mean((fhat - indicator)^2)
        }
        brier[j] <- sum(brier_weights * bs_c)
      }
      list(brier = brier)
    } else if (type == "loglik") {
      mass <- predict_loglik_mass_cpmnet(fit_fold, newx = x_val,
                                         y_test = y_val, s = NULL)
      negll <- rep(NA_real_, nlambda)
      for (j in 1:nlambda_fold) {
        negll[j] <- -mean(log(pmax(mass[, j], 1e-12)))
      }
      list(negll = negll)
    } else {
      # psr_abs or psr_sq
      cdf_vals <- predict_psr_cpmnet(fit_fold, newx = x_val,
                                      y_test = y_val, s = NULL)
      psr <- rep(NA_real_, nlambda)
      for (j in 1:nlambda_fold) {
        psr_vals <- 2.0 * cdf_vals[, j] - 1.0
        if (type == "psr_abs") {
          psr[j] <- mean(abs(psr_vals))
        } else {
          psr[j] <- mean(psr_vals^2)
        }
      }
      list(psr = psr)
    }
  }
  max_fold <- max(foldid)
  if (parallel && max_fold > 1 && .Platform$OS.type != "windows") {
    if (is.null(ncores)) ncores <- min(max_fold, parallel::detectCores() - 1)
    fold_results <- parallel::mclapply(seq(max_fold), process_fold, mc.cores = ncores)
  } else {
    if (parallel && .Platform$OS.type == "windows") {
      warning("Parallel processing not supported on Windows. Running sequentially.",
              call. = FALSE)
    }
    fold_results <- lapply(seq(max_fold), process_fold)
  }
  has_error <- sapply(fold_results, inherits, "try-error")
  if (any(has_error)) {
    n_failed <- sum(has_error)
    warning(paste(n_failed, "fold(s) failed during parallel execution."),
            call. = FALSE)
    fold_results <- fold_results[!has_error]
  }
  if (length(fold_results) == 0) stop("All folds failed.", call. = FALSE)
  nfolds_actual <- length(fold_results)
  cv_call <- match.call()
  if (type %in% c("mean", "median")) {
    mae <- .cv_min_1se(do.call(rbind, lapply(fold_results, `[[`, "mae")))
    mse <- .cv_min_1se(do.call(rbind, lapply(fold_results, `[[`, "mse")))
    pR2 <- .cv_min_1se(do.call(rbind, lapply(fold_results, `[[`, "pR2")),
                       maximize = TRUE, allow_empty = TRUE)
    nzero <- fit_full$df
    names(nzero) <- paste0("s", round(lambda_path, 7))
    index <- matrix(c(mae$opt_idx, mae$ise_idx, mse$opt_idx, mse$ise_idx,
                      pR2$opt_idx, pR2$ise_idx),
                    nrow = 2, ncol = 3,
                    dimnames = list(c("min", "1se"), c("MAE", "MSE", "pR2")))
    structure(
      list(
        lambda = lambda_path,
        cvm.mae = mae$mean, cvsd.mae = mae$se,
        cvup.mae = mae$mean + mae$se, cvlo.mae = mae$mean - mae$se,
        lambda.min.mae = lambda_path[mae$opt_idx],
        lambda.1se.mae = lambda_path[mae$ise_idx],
        cvm.mse = mse$mean, cvsd.mse = mse$se,
        cvup.mse = mse$mean + mse$se, cvlo.mse = mse$mean - mse$se,
        lambda.min.mse = lambda_path[mse$opt_idx],
        lambda.1se.mse = lambda_path[mse$ise_idx],
        cvm.pR2 = pR2$mean, cvsd.pR2 = pR2$se,
        cvup.pR2 = pR2$mean + pR2$se, cvlo.pR2 = pR2$mean - pR2$se,
        lambda.max.pR2 = if (!is.na(pR2$opt_idx)) lambda_path[pR2$opt_idx] else NA,
        lambda.1se.pR2 = if (!is.na(pR2$ise_idx)) lambda_path[pR2$ise_idx] else NA,
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type, nfolds_successful = nfolds_actual, call = cv_call
      ),
      class = "cv.cpmnet"
    )
  } else {
    # single criterion types
    slot <- switch(type,
                   pinball_abs = "pinball", pinball_sq = "pinball",
                   brier = "brier", loglik = "negll", "psr")
    extra <- switch(type,
                    pinball_abs = ,
                    pinball_sq = list(tau_levels = tau_levels,
                                      tau_weights = tau_weights),
                    brier = list(brier_probs = brier_probs,
                                 brier_weights = brier_weights,
                                 brier_cuts = brier_cuts),
                    NULL)
    sel <- .cv_min_1se(do.call(rbind, lapply(fold_results, `[[`, slot)))
    .cv_result(lambda_path, sel, fit_full, type, extra, nfolds_actual, cv_call)
  }
}

# validate a probability grid and weights
#' @keywords internal
#' @noRd
.check_prob_grid <- function(values, weights, val_name, wt_name) {
  if (!is.numeric(values) || length(values) == 0L) {
    stop(val_name, " must be a non-empty numeric vector; got ",
         deparse(values), call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop(val_name, " must contain no NA/NaN/Inf", call. = FALSE)
  }
  if (any(values <= 0) || any(values >= 1)) {
    stop(val_name, " must lie strictly in (0, 1); got range [",
         min(values), ", ", max(values), "]", call. = FALSE)
  }
  if (is.unsorted(values)) {
    stop(val_name, " must be sorted in increasing order", call. = FALSE)
  }
  n_val <- length(values)
  if (is.null(weights)) {
    weights <- rep(1.0 / n_val, n_val)
  } else {
    if (!is.numeric(weights) || length(weights) != n_val) {
      stop(wt_name, " must be numeric of length ", n_val,
           " (matching ", val_name, "); got length ", length(weights),
           call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop(wt_name, " must be non-negative and finite", call. = FALSE)
    }
    if (abs(sum(weights) - 1.0) > 1e-8) {
      stop(wt_name, " must sum to 1; got sum = ", sum(weights),
           call. = FALSE)
    }
  }
  weights
}

# min and 1se lambda indices
#' @keywords internal
#' @noRd
.cv_min_1se <- function(fold_mat, maximize = FALSE, allow_empty = FALSE) {
  m <- colMeans(fold_mat, na.rm = TRUE)
  se <- apply(fold_mat, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_mat)))
  valid <- which(is.finite(m))
  if (length(valid) == 0) {
    if (allow_empty) return(list(mean = m, se = se, opt_idx = NA, ise_idx = NA))
    stop("No valid CV results.", call. = FALSE)
  }
  opt_idx <- if (maximize) valid[which.max(m[valid])]
             else valid[which.min(m[valid])]
  if (maximize) {
    threshold <- m[opt_idx] - se[opt_idx]
    cand <- which(m >= threshold & is.finite(m))
  } else {
    threshold <- m[opt_idx] + se[opt_idx]
    cand <- which(m <= threshold & is.finite(m))
  }
  ise_idx <- if (length(cand) == 0) opt_idx else min(cand)
  list(mean = m, se = se, opt_idx = opt_idx, ise_idx = ise_idx)
}

# assemble the return object
#' @keywords internal
#' @noRd
.cv_result <- function(lambda_path, sel, fit_full, type, extra,
                       nfolds_actual, call) {
  nzero <- fit_full$df
  names(nzero) <- paste0("s", round(lambda_path, 7))
  index <- matrix(c(sel$opt_idx, sel$ise_idx), nrow = 2, ncol = 1,
                  dimnames = list(c("min", "1se"), type))
  structure(
    c(list(
        lambda = lambda_path,
        cvm = sel$mean, cvsd = sel$se,
        cvup = sel$mean + sel$se, cvlo = sel$mean - sel$se,
        lambda.min = lambda_path[sel$opt_idx],
        lambda.1se = lambda_path[sel$ise_idx],
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type),
      extra,
      list(nfolds_successful = nfolds_actual, call = call)),
    class = "cv.cpmnet"
  )
}
