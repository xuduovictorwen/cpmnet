#' Plot coefficient paths from a cpmnet fit
#'
#' Produces a glmnet-style coefficient path plot. Each line traces one
#' predictor's slope across the penalization (regularization) path. A top-axis shows
#' the number of nonzero coefficients.
#'
#' @param x fitted \code{"cpmnet"} object.
#' @param xvar x-axis variable. One of \code{"lambda"} (log-lambda),
#'   \code{"norm"} (L1 norm of the slope vector), or \code{"dev"}
#'   (fraction deviance explained).
#' @param label if \code{TRUE}, label nonzero coefficients at the
#'   rightmost point with their predictor index.
#' @param ... passed to \code{\link[graphics]{plot}}.
#'
#' @return \code{NULL}, invisibly. Called for its side effect.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{plot.cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' fit <- cpmnet(x, y, family = "probit", nlambda = 20)
#' plot(fit, xvar = "lambda")
#' }
#'
#' @export
plot.cpmnet <- function(x, xvar = c("lambda", "norm", "dev"),
                        label = FALSE, ...) {

  xvar <- match.arg(xvar)
  beta <- x$beta
  lambda <- x$lambda
  df <- x$df
  p <- nrow(beta)

  xval <- switch(xvar,
                 lambda = log(lambda),
                 norm = colSums(abs(beta)),
                 dev = x$dev.ratio)

  xlab <- switch(xvar,
                 lambda = "Log Lambda",
                 norm = "L1 Norm",
                 dev = "Fraction Deviance Explained")

  ylim <- range(beta)
  plot(xval, beta[1, ], type = "n", xlim = rev(range(xval)), ylim = ylim,
       xlab = xlab, ylab = "Coefficients", ...)

  matlines(xval, t(beta), lty = 1, lwd = 1.5, col = rainbow(p))

  if (label) {
    nonzero <- which(abs(beta[, ncol(beta)]) > 0)
    if (length(nonzero) > 0)
      text(max(xval), beta[nonzero, ncol(beta)],
           labels = nonzero, cex = 0.7, pos = 4)
  }

  n_ticks <- min(10, length(xval))
  tick_idx <- seq(1, length(xval), length.out = n_ticks)
  axis(3, at = xval[tick_idx], labels = df[tick_idx], cex.axis = 0.8)
  mtext("Degrees of Freedom", side = 3, line = 2.5, cex = 0.8)
  grid()
  invisible(NULL)
}

#' Plot cross-validation curves from cv.cpmnet
#'
#' Plots the selected CV metric (mean across folds) against log-lambda,
#' with vertical bars showing plus/minus one standard error. Dashed
#' lines mark \code{lambda.min} (black) and \code{lambda.1se} (blue).
#'
#' @param x fitted \code{"cv.cpmnet"} object.
#' @param metric which metric to plot. Only used when \code{x$type} is
#'   \code{"mean"} or \code{"median"}; one of \code{"mae"} (default),
#'   \code{"mse"}, or \code{"pR2"}. Ignored for other types.
#' @param sign.lambda \code{1} to plot against \code{log(lambda)} (the
#'   default), \code{-1} to flip the x-axis.
#' @param ... passed to \code{\link[graphics]{plot}}.
#'
#' @return \code{NULL}, invisibly.
#'
#' @seealso \code{\link{cv.cpmnet}}, \code{\link{plot.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(120 * 5), 120, 5)
#' y <- x[, 1] - x[, 2] + rnorm(120)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5, type = "brier",
#'                    nlambda = 15, parallel = FALSE)
#' plot(cvfit)
#' }
#'
#' @export
plot.cv.cpmnet <- function(x, metric = NULL, sign.lambda = 1, ...) {
  if (x$type %in% c("pinball_abs", "pinball_sq", "brier", "psr_abs", "psr_sq",
                    "loglik")) {
    cvm <- x$cvm
    cvsd <- x$cvsd
    cvup <- x$cvup
    cvlo <- x$cvlo
    lambda_opt <- x$lambda.min
    lambda_1se <- x$lambda.1se
    ylab <- switch(x$type,
                   pinball_abs = "Pinball Loss (Absolute)",
                   pinball_sq = "Pinball Loss (Squared)",
                   brier = "Brier Score",
                   psr_abs = "PSR (Absolute)",
                   psr_sq = "PSR (Squared)",
                   loglik = "Negative Log-Likelihood")
  } else {
    if (is.null(metric)) metric <- "mae"
    valid_metrics <- c("mae", "mse", "pR2")
    if (!metric %in% valid_metrics) {
      stop("metric for type '", x$type, "' must be one of: ",
           paste(valid_metrics, collapse = ", "), "; got '", metric, "'",
           call. = FALSE)
    }
    cvm <- x[[paste0("cvm.", metric)]]
    cvsd <- x[[paste0("cvsd.", metric)]]
    cvup <- x[[paste0("cvup.", metric)]]
    cvlo <- x[[paste0("cvlo.", metric)]]
    if (metric == "pR2") {
      lambda_opt <- x$lambda.max.pR2
      lambda_1se <- x$lambda.1se.pR2
    } else {
      lambda_opt <- x[[paste0("lambda.min.", metric)]]
      lambda_1se <- x[[paste0("lambda.1se.", metric)]]
    }
    ylab <- switch(metric,
                   mae = "Mean Absolute Error",
                   mse = "Mean-Squared Error",
                   pR2 = "pR2 (PA precision * PA accuracy)")
  }
  lambda <- x$lambda
  nzero <- x$nzero
  xval <- sign.lambda * log(lambda)
  xlab <- if (sign.lambda < 0) "-Log(Lambda)" else "Log(Lambda)"
  ylim <- range(cvlo, cvup, na.rm = TRUE)
  plot(xval, cvm, type = "n", xlim = rev(range(xval)), ylim = ylim,
       xlab = xlab, ylab = ylab, ...)
  finite_idx <- which(is.finite(cvm) & is.finite(cvlo) & is.finite(cvup))
  if (length(finite_idx) > 0) {
    arrows(xval[finite_idx], cvlo[finite_idx], xval[finite_idx], cvup[finite_idx],
           angle = 90, code = 3, length = 0.05, col = "black", lwd = 0.5)
    points(xval[finite_idx], cvm[finite_idx], pch = 20, col = "red", cex = 0.8)
  }
  if (!is.na(lambda_opt))
    abline(v = sign.lambda * log(lambda_opt), lty = 2, col = "black", lwd = 0.5)
  if (!is.na(lambda_1se) && abs(lambda_opt - lambda_1se) > 0)
    abline(v = sign.lambda * log(lambda_1se), lty = 2, col = "blue", lwd = 0.5)
  n_ticks <- min(10, length(xval))
  tick_idx <- seq(1, length(xval), length.out = n_ticks)
  axis(3, at = xval[tick_idx], labels = nzero[tick_idx], cex.axis = 0.8)
  mtext("Degrees of Freedom", side = 3, line = 2.5, cex = 0.8)
  legend("topright", legend = c("min", "1se"), col = c("black", "blue"),
         lty = 2, lwd = 0.8, cex = 0.7, bty = "o", bg = "white")
  grid()
  invisible(NULL)
}
