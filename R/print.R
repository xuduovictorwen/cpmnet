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
    cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
    cat(sprintf("%5s %7.5f %5d %7.3f %7.5f %7d\n",
                "min", x$lambda.min.mae, x$index["min", "MAE"],
                x$cvm.mae[x$index["min", "MAE"]], x$cvsd.mae[x$index["min", "MAE"]],
                x$nzero[x$index["min", "MAE"]]))
    cat(sprintf("%5s %7.5f %5d %7.3f %7.5f %7d\n",
                "1se", x$lambda.1se.mae, x$index["1se", "MAE"],
                x$cvm.mae[x$index["1se", "MAE"]], x$cvsd.mae[x$index["1se", "MAE"]],
                x$nzero[x$index["1se", "MAE"]]))
    cat("\nMeasure: Mean-Squared Error \n\n")
    cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
    cat(sprintf("%5s %7.5f %5d %7.3f %7.5f %7d\n",
                "min", x$lambda.min.mse, x$index["min", "MSE"],
                x$cvm.mse[x$index["min", "MSE"]], x$cvsd.mse[x$index["min", "MSE"]],
                x$nzero[x$index["min", "MSE"]]))
    cat(sprintf("%5s %7.5f %5d %7.3f %7.5f %7d\n",
                "1se", x$lambda.1se.mse, x$index["1se", "MSE"],
                x$cvm.mse[x$index["1se", "MSE"]], x$cvsd.mse[x$index["1se", "MSE"]],
                x$nzero[x$index["1se", "MSE"]]))
    if (x$type == "mean" && !is.na(x$lambda.max.pR2)) {
      cat("\nMeasure: pR2 (PA precision * PA accuracy) \n\n")
      cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
      cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                  "max", x$lambda.max.pR2, x$index["min", "pR2"],
                  x$cvm.pR2[x$index["min", "pR2"]], x$cvsd.pR2[x$index["min", "pR2"]],
                  x$nzero[x$index["min", "pR2"]]))
      cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                  "1se", x$lambda.1se.pR2, x$index["1se", "pR2"],
                  x$cvm.pR2[x$index["1se", "pR2"]], x$cvsd.pR2[x$index["1se", "pR2"]],
                  x$nzero[x$index["1se", "pR2"]]))
    }
  } else if (x$type %in% c("pinball_abs", "pinball_sq")) {
    label <- if (x$type == "pinball_abs") "Pinball Loss (Absolute)" else "Pinball Loss (Squared)"
    cat(sprintf("Measure: %s\n", label))
    cat(sprintf("  tau levels: %s\n\n", paste(x$tau_levels, collapse = ", ")))
    cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "min", x$lambda.min, x$index["min", 1],
                x$cvm[x$index["min", 1]], x$cvsd[x$index["min", 1]],
                x$nzero[x$index["min", 1]]))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "1se", x$lambda.1se, x$index["1se", 1],
                x$cvm[x$index["1se", 1]], x$cvsd[x$index["1se", 1]],
                x$nzero[x$index["1se", 1]]))
  } else if (x$type %in% c("psr_abs", "psr_sq")) {
    label <- if (x$type == "psr_abs") "PSR (Absolute)" else "PSR (Squared)"
    cat(sprintf("Measure: %s\n\n", label))
    cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "min", x$lambda.min, x$index["min", 1],
                x$cvm[x$index["min", 1]], x$cvsd[x$index["min", 1]],
                x$nzero[x$index["min", 1]]))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "1se", x$lambda.1se, x$index["1se", 1],
                x$cvm[x$index["1se", 1]], x$cvsd[x$index["1se", 1]],
                x$nzero[x$index["1se", 1]]))
  } else {
    cat("Measure: Brier Score\n")
    cat(sprintf("  thresholds (percentiles): %s\n", paste(x$brier_probs, collapse = ", ")))
    cat(sprintf("  threshold values: %s\n\n", paste(round(x$brier_cuts, 4), collapse = ", ")))
    cat(sprintf("%5s %7s %5s %7s %7s %7s\n", "", "Lambda", "Index", "Measure", "SE", "Nonzero"))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "min", x$lambda.min, x$index["min", 1],
                x$cvm[x$index["min", 1]], x$cvsd[x$index["min", 1]],
                x$nzero[x$index["min", 1]]))
    cat(sprintf("%5s %7.5f %5d %7.4f %7.5f %7d\n",
                "1se", x$lambda.1se, x$index["1se", 1],
                x$cvm[x$index["1se", 1]], x$cvsd[x$index["1se", 1]],
                x$nzero[x$index["1se", 1]]))
  }
  cat("\n")
  invisible(x)
}
