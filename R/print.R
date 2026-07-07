#' Print a cpmnet fit
#'
#' Prints a glmnet-style summary of a \code{cpmnet} fit: the function
#' call followed by a table with per-lambda degrees of freedom,
#' percent deviance explained, and lambda value.
#'
#' @param x fitted \code{"cpmnet"} object.
#' @param ... ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{summary.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' fit <- cpmnet(x, y, family = "probit", nlambda = 10)
#' print(fit)
#' }
#'
#' @export
print.cpmnet <- function(x, ...) {
  family_name <- names(rms::probabilityFamilies)[x$data_prep$link]
  cat("Call:  cpmnet(x = x, y = y, family = \"", family_name,
      "\", alpha_en = ", x$alpha_en,
      ", nlambda = ", length(x$lambda), ")\n\n", sep = "")

  cat(sprintf("%4s %3s %6s %10s\n", "", "Df", "%Dev", "Lambda"))
  for (i in seq_along(x$lambda)) {
    cat(sprintf("%3d %3d %6.2f %20.15f\n",
                i, x$df[i], x$dev.ratio[i] * 100, x$lambda[i]))
  }
  invisible(x)
}

#' Print a cv.cpmnet fit
#'
#' Prints the CV metric selected by \code{type} at \code{lambda.min} and
#' \code{lambda.1se}. For \code{type} in \code{"mean"}/\code{"median"},
#' MAE and MSE (and pR2 when \code{type = "mean"}) are reported.
#'
#' @param x fitted \code{"cv.cpmnet"} object.
#' @param ... ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{cv.cpmnet}}, \code{\link{summary.cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5, type = "brier",
#'                    nlambda = 10, parallel = FALSE)
#' print(cvfit)
#' }
#'
#' @export
print.cv.cpmnet <- function(x, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  if (x$type %in% c("mean", "median")) {
    cat("Measure: Mean Absolute Error \n\n")
    .print_cv_rows(c("min", "1se"),
                   c(x$lambda.min.mae, x$lambda.1se.mae),
                   x$index[, "MAE"], x$cvm.mae, x$cvsd.mae, x$nzero, "%7.3f")
    cat("\nMeasure: Mean-Squared Error \n\n")
    .print_cv_rows(c("min", "1se"),
                   c(x$lambda.min.mse, x$lambda.1se.mse),
                   x$index[, "MSE"], x$cvm.mse, x$cvsd.mse, x$nzero, "%7.3f")
    if (x$type == "mean" && !is.na(x$lambda.max.pR2)) {
      cat("\nMeasure: pR2 (PA precision * PA accuracy) \n\n")
      .print_cv_rows(c("max", "1se"),
                     c(x$lambda.max.pR2, x$lambda.1se.pR2),
                     x$index[, "pR2"], x$cvm.pR2, x$cvsd.pR2, x$nzero)
    }
  } else {
    header <- switch(x$type,
      pinball_abs = "Measure: Pinball Loss (Absolute)\n",
      pinball_sq  = "Measure: Pinball Loss (Squared)\n",
      loglik      = "Measure: Negative Log-Likelihood (held-out)\n\n",
      psr_abs     = "Measure: PSR (Absolute)\n\n",
      psr_sq      = "Measure: PSR (Squared)\n\n",
      brier       = "Measure: Brier Score\n")
    cat(header)
    if (x$type %in% c("pinball_abs", "pinball_sq")) {
      cat(sprintf("  tau levels: %s\n\n", paste(x$tau_levels, collapse = ", ")))
    } else if (x$type == "brier") {
      cat(sprintf("  thresholds (percentiles): %s\n", paste(x$brier_probs, collapse = ", ")))
      cat(sprintf("  threshold values: %s\n\n", paste(round(x$brier_cuts, 4), collapse = ", ")))
    }
    .print_cv_rows(c("min", "1se"), c(x$lambda.min, x$lambda.1se),
                   x$index[, 1], x$cvm, x$cvsd, x$nzero)
  }
  cat("\n")
  invisible(x)
}

# Two-row lambda-selection table shared by the print.cv.cpmnet branches.
# idx holds the (min/opt, 1se) lambda indices; measure_fmt controls the
# Measure column format (mean/median print %7.3f).
#' @keywords internal
#' @noRd
.print_cv_rows <- function(labels, lambdas, idx, cvm, cvsd, nzero,
                           measure_fmt = "%7.4f") {
  cat(sprintf("%5s %7s %5s %7s %7s %7s\n",
              "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
  row_fmt <- paste0("%5s %7.5f %5d ", measure_fmt, " %7.5f %7d\n")
  for (r in 1:2) {
    cat(sprintf(row_fmt, labels[r], lambdas[r], idx[r],
                cvm[idx[r]], cvsd[idx[r]], nzero[idx[r]]))
  }
}
