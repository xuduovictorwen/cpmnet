#' Cross-validation for cpmnet
#'
#' Performs k-fold cross-validation for a \code{cpmnet} model and selects
#' lambda according to one of four conceptual criteria: mean, median,
#' pinball, or Brier. The default, \code{"brier"}, averages the squared
#' Brier score across a user-specified grid of response percentiles and
#' is the criterion recommended in the accompanying methods paper.
#'
#' @param x predictor matrix. See \code{\link{cpmnet}} for requirements.
#' @param y response vector. See \code{\link{cpmnet}}.
#' @param nfolds number of folds. Must satisfy \code{2 <= nfolds <= n}.
#' @param foldid optional integer vector of length \code{n} assigning each
#'   observation to a fold. When supplied, values must form
#'   \code{1:max(foldid)} with no gaps.
#' @param type cross-validation loss. One of \code{"brier"} (default),
#'   \code{"mean"}, \code{"median"}, \code{"pinball_abs"},
#'   \code{"pinball_sq"}, \code{"psr_abs"}, \code{"psr_sq"}.
#'   These map to four conceptual criteria emphasized in the paper:
#'   mean and median regression loss (with MAE/MSE/pR2 reported),
#'   pinball quantile loss (absolute or squared variant), and the
#'   Brier score. The \code{"psr_*"} options remain in the package
#'   for completeness but are not part of the paper's primary set.
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
                               "pinball_sq", "psr_abs", "psr_sq"),
                      parallel = TRUE, ncores = NULL,
                      tau_levels = seq(0.1, 0.9, by = 0.1),
                      tau_weights = NULL,
                      brier_probs = seq(0.1, 0.9, by = 0.1),
                      brier_weights = NULL, ...) {

  type <- match.arg(type)

  # dimensional baseline (so validation errors can reference n)
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

  # tau_levels (used by pinball types)
  if (!is.numeric(tau_levels) || length(tau_levels) == 0L) {
    stop("tau_levels must be a non-empty numeric vector; got ",
         deparse(tau_levels), call. = FALSE)
  }
  if (any(!is.finite(tau_levels))) {
    stop("tau_levels must contain no NA/NaN/Inf", call. = FALSE)
  }
  if (any(tau_levels <= 0) || any(tau_levels >= 1)) {
    stop("tau_levels must lie strictly in (0, 1); got range [",
         min(tau_levels), ", ", max(tau_levels), "]", call. = FALSE)
  }
  if (is.unsorted(tau_levels)) {
    stop("tau_levels must be sorted in increasing order", call. = FALSE)
  }
  n_tau <- length(tau_levels)
  if (is.null(tau_weights)) {
    tau_weights <- rep(1.0 / n_tau, n_tau)
  } else {
    if (!is.numeric(tau_weights) || length(tau_weights) != n_tau) {
      stop("tau_weights must be numeric of length ", n_tau,
           " (matching tau_levels); got length ", length(tau_weights),
           call. = FALSE)
    }
    if (any(!is.finite(tau_weights)) || any(tau_weights < 0)) {
      stop("tau_weights must be non-negative and finite", call. = FALSE)
    }
    if (abs(sum(tau_weights) - 1.0) > 1e-8) {
      stop("tau_weights must sum to 1; got sum = ", sum(tau_weights),
           call. = FALSE)
    }
  }

  # brier_probs (used by brier type)
  if (!is.numeric(brier_probs) || length(brier_probs) == 0L) {
    stop("brier_probs must be a non-empty numeric vector; got ",
         deparse(brier_probs), call. = FALSE)
  }
  if (any(!is.finite(brier_probs))) {
    stop("brier_probs must contain no NA/NaN/Inf", call. = FALSE)
  }
  if (any(brier_probs <= 0) || any(brier_probs >= 1)) {
    stop("brier_probs must lie strictly in (0, 1); got range [",
         min(brier_probs), ", ", max(brier_probs), "]", call. = FALSE)
  }
  if (is.unsorted(brier_probs)) {
    stop("brier_probs must be sorted in increasing order", call. = FALSE)
  }
  n_brier <- length(brier_probs)
  if (is.null(brier_weights)) {
    brier_weights <- rep(1.0 / n_brier, n_brier)
  } else {
    if (!is.numeric(brier_weights) || length(brier_weights) != n_brier) {
      stop("brier_weights must be numeric of length ", n_brier,
           " (matching brier_probs); got length ", length(brier_weights),
           call. = FALSE)
    }
    if (any(!is.finite(brier_weights)) || any(brier_weights < 0)) {
      stop("brier_weights must be non-negative and finite", call. = FALSE)
    }
    if (abs(sum(brier_weights) - 1.0) > 1e-8) {
      stop("brier_weights must sum to 1; got sum = ", sum(brier_weights),
           call. = FALSE)
    }
  }

  dots <- list(...)
  fit_full <- do.call(cpmnet, c(list(x = x, y = y), dots))
  lambda_path <- fit_full$lambda
  nlambda <- length(lambda_path)
  link_name <- names(rms::probabilityFamilies)[fit_full$data_prep$link]
  cumprob <- eval(rms::probabilityFamilies[[link_name]][1])
  # brier thresholds from all y before splitting
  brier_cuts <- as.numeric(quantile(y, probs = brier_probs))

  process_fold <- function(fold) {
    val_idx <- which(foldid == fold)
    train_idx <- which(foldid != fold)
    n_val <- length(val_idx)
    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_val <- x[val_idx, , drop = FALSE]
    y_val <- y[val_idx]
    fit_fold <- do.call(cpmnet, c(list(x = x_train, y = y_train,
                                       lambda = lambda_path,
                                       verbose = FALSE), dots))
    nlambda_fold <- length(fit_fold$lambda)
    if (type %in% c("mean", "median")) {
      # predict.cpmnet always returns a matrix n_val x nlambda_fold
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
  if (type %in% c("mean", "median")) {
    fold_mae <- do.call(rbind, lapply(fold_results, `[[`, "mae"))
    fold_mse <- do.call(rbind, lapply(fold_results, `[[`, "mse"))
    fold_pR2 <- do.call(rbind, lapply(fold_results, `[[`, "pR2"))
    mae_mean <- colMeans(fold_mae, na.rm = TRUE)
    mae_se <- apply(fold_mae, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_mae)))
    mse_mean <- colMeans(fold_mse, na.rm = TRUE)
    mse_se <- apply(fold_mse, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_mse)))
    pR2_mean <- colMeans(fold_pR2, na.rm = TRUE)
    pR2_se <- apply(fold_pR2, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_pR2)))
    valid_mae <- which(is.finite(mae_mean))
    if (length(valid_mae) == 0) stop("No valid CV results.", call. = FALSE)
    mae_min_idx <- valid_mae[which.min(mae_mean[valid_mae])]
    valid_mse <- which(is.finite(mse_mean))
    mse_min_idx <- valid_mse[which.min(mse_mean[valid_mse])]
    valid_pR2 <- which(is.finite(pR2_mean))
    pR2_max_idx <- if (length(valid_pR2) > 0) valid_pR2[which.max(pR2_mean[valid_pR2])] else NA
    mae_threshold <- mae_mean[mae_min_idx] + mae_se[mae_min_idx]
    mae_1se_candidates <- which(mae_mean <= mae_threshold & is.finite(mae_mean))
    mae_1se_idx <- if (length(mae_1se_candidates) == 0) mae_min_idx else min(mae_1se_candidates)
    mse_threshold <- mse_mean[mse_min_idx] + mse_se[mse_min_idx]
    mse_1se_candidates <- which(mse_mean <= mse_threshold & is.finite(mse_mean))
    mse_1se_idx <- if (length(mse_1se_candidates) == 0) mse_min_idx else min(mse_1se_candidates)
    if (!is.na(pR2_max_idx)) {
      pR2_threshold <- pR2_mean[pR2_max_idx] - pR2_se[pR2_max_idx]
      pR2_1se_candidates <- which(pR2_mean >= pR2_threshold & is.finite(pR2_mean))
      pR2_1se_idx <- if (length(pR2_1se_candidates) == 0) pR2_max_idx else min(pR2_1se_candidates)
    } else {
      pR2_1se_idx <- NA
    }
    nzero <- fit_full$df
    names(nzero) <- paste0("s", round(lambda_path, 7))
    index <- matrix(c(mae_min_idx, mae_1se_idx, mse_min_idx, mse_1se_idx,
                       pR2_max_idx, pR2_1se_idx),
                    nrow = 2, ncol = 3,
                    dimnames = list(c("min", "1se"), c("MAE", "MSE", "pR2")))
    structure(
      list(
        lambda = lambda_path,
        cvm.mae = mae_mean, cvsd.mae = mae_se,
        cvup.mae = mae_mean + mae_se, cvlo.mae = mae_mean - mae_se,
        lambda.min.mae = lambda_path[mae_min_idx],
        lambda.1se.mae = lambda_path[mae_1se_idx],
        cvm.mse = mse_mean, cvsd.mse = mse_se,
        cvup.mse = mse_mean + mse_se, cvlo.mse = mse_mean - mse_se,
        lambda.min.mse = lambda_path[mse_min_idx],
        lambda.1se.mse = lambda_path[mse_1se_idx],
        cvm.pR2 = pR2_mean, cvsd.pR2 = pR2_se,
        cvup.pR2 = pR2_mean + pR2_se, cvlo.pR2 = pR2_mean - pR2_se,
        lambda.max.pR2 = if (!is.na(pR2_max_idx)) lambda_path[pR2_max_idx] else NA,
        lambda.1se.pR2 = if (!is.na(pR2_1se_idx)) lambda_path[pR2_1se_idx] else NA,
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type, nfolds_successful = nfolds_actual, call = match.call()
      ),
      class = "cv.cpmnet"
    )
  } else if (type %in% c("pinball_abs", "pinball_sq")) {
    fold_pb <- do.call(rbind, lapply(fold_results, `[[`, "pinball"))
    pb_mean <- colMeans(fold_pb, na.rm = TRUE)
    pb_se <- apply(fold_pb, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_pb)))
    valid_pb <- which(is.finite(pb_mean))
    if (length(valid_pb) == 0) stop("No valid CV results.", call. = FALSE)
    pb_min_idx <- valid_pb[which.min(pb_mean[valid_pb])]
    pb_threshold <- pb_mean[pb_min_idx] + pb_se[pb_min_idx]
    pb_1se_candidates <- which(pb_mean <= pb_threshold & is.finite(pb_mean))
    pb_1se_idx <- if (length(pb_1se_candidates) == 0) pb_min_idx else min(pb_1se_candidates)
    nzero <- fit_full$df
    names(nzero) <- paste0("s", round(lambda_path, 7))
    index <- matrix(c(pb_min_idx, pb_1se_idx), nrow = 2, ncol = 1,
                    dimnames = list(c("min", "1se"), type))
    structure(
      list(
        lambda = lambda_path,
        cvm = pb_mean, cvsd = pb_se,
        cvup = pb_mean + pb_se, cvlo = pb_mean - pb_se,
        lambda.min = lambda_path[pb_min_idx],
        lambda.1se = lambda_path[pb_1se_idx],
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type, tau_levels = tau_levels, tau_weights = tau_weights,
        nfolds_successful = nfolds_actual, call = match.call()
      ),
      class = "cv.cpmnet"
    )
  } else if (type == "brier") {
    fold_bs <- do.call(rbind, lapply(fold_results, `[[`, "brier"))
    bs_mean <- colMeans(fold_bs, na.rm = TRUE)
    bs_se <- apply(fold_bs, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_bs)))
    valid_bs <- which(is.finite(bs_mean))
    if (length(valid_bs) == 0) stop("No valid CV results.", call. = FALSE)
    bs_min_idx <- valid_bs[which.min(bs_mean[valid_bs])]
    bs_threshold <- bs_mean[bs_min_idx] + bs_se[bs_min_idx]
    bs_1se_candidates <- which(bs_mean <= bs_threshold & is.finite(bs_mean))
    bs_1se_idx <- if (length(bs_1se_candidates) == 0) bs_min_idx else min(bs_1se_candidates)
    nzero <- fit_full$df
    names(nzero) <- paste0("s", round(lambda_path, 7))
    index <- matrix(c(bs_min_idx, bs_1se_idx), nrow = 2, ncol = 1,
                    dimnames = list(c("min", "1se"), "brier"))
    structure(
      list(
        lambda = lambda_path,
        cvm = bs_mean, cvsd = bs_se,
        cvup = bs_mean + bs_se, cvlo = bs_mean - bs_se,
        lambda.min = lambda_path[bs_min_idx],
        lambda.1se = lambda_path[bs_1se_idx],
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type, brier_probs = brier_probs, brier_weights = brier_weights,
        brier_cuts = brier_cuts,
        nfolds_successful = nfolds_actual, call = match.call()
      ),
      class = "cv.cpmnet"
    )
  } else {
    # psr_abs or psr_sq
    fold_psr <- do.call(rbind, lapply(fold_results, `[[`, "psr"))
    psr_mean <- colMeans(fold_psr, na.rm = TRUE)
    psr_se <- apply(fold_psr, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(fold_psr)))
    valid_psr <- which(is.finite(psr_mean))
    if (length(valid_psr) == 0) stop("No valid CV results.", call. = FALSE)
    psr_min_idx <- valid_psr[which.min(psr_mean[valid_psr])]
    psr_threshold <- psr_mean[psr_min_idx] + psr_se[psr_min_idx]
    psr_1se_candidates <- which(psr_mean <= psr_threshold & is.finite(psr_mean))
    psr_1se_idx <- if (length(psr_1se_candidates) == 0) psr_min_idx else min(psr_1se_candidates)
    nzero <- fit_full$df
    names(nzero) <- paste0("s", round(lambda_path, 7))
    index <- matrix(c(psr_min_idx, psr_1se_idx), nrow = 2, ncol = 1,
                    dimnames = list(c("min", "1se"), type))
    structure(
      list(
        lambda = lambda_path,
        cvm = psr_mean, cvsd = psr_se,
        cvup = psr_mean + psr_se, cvlo = psr_mean - psr_se,
        lambda.min = lambda_path[psr_min_idx],
        lambda.1se = lambda_path[psr_1se_idx],
        nzero = nzero, cpmnet.fit = fit_full, index = index,
        type = type,
        nfolds_successful = nfolds_actual, call = match.call()
      ),
      class = "cv.cpmnet"
    )
  }
}
