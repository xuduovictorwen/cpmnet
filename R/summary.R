#' Summarize a cpmnet fit
#'
#' Tabulates lambda, df, fraction of deviance explained, and convergence
#' flag for each lambda in the fitted path.
#'
#' @param object fitted \code{"cpmnet"} object.
#' @param ... ignored.
#'
#' @return A data frame with one row per lambda, invisibly.
#'
#' @seealso \code{\link{cpmnet}}, \code{\link{print.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' fit <- cpmnet(x, y, family = "probit", nlambda = 10)
#' summary(fit)
#' }
#'
#' @export
summary.cpmnet <- function(object, ...) {
  if (!inherits(object, "cpmnet")) {
    stop("object must inherit from class 'cpmnet'", call. = FALSE)
  }
  tbl <- data.frame(
    index      = seq_along(object$lambda),
    lambda     = object$lambda,
    df         = object$df,
    dev.ratio  = object$dev.ratio,
    logL       = object$logL,
    converged  = object$converged,
    niter      = object$niter
  )
  cat("cpmnet fit:\n")
  cat(sprintf("  nobs        : %d\n", object$nobs))
  cat(sprintf("  nlambda     : %d\n", length(object$lambda)))
  cat(sprintf("  alpha_en    : %.3f\n", object$alpha_en))
  cat(sprintf("  all converged: %s\n",
              if (all(object$converged)) "yes" else "no"))
  cat(sprintf("  max df      : %d\n", max(object$df)))
  cat(sprintf("  max dev.ratio: %.4f\n", max(object$dev.ratio)))
  cat("\nPath (head):\n")
  print(utils::head(tbl, 10))
  invisible(tbl)
}

#' Summarize a cv.cpmnet fit
#'
#' Reports the CV criterion, the chosen \code{lambda.min} and
#' \code{lambda.1se}, and the cv metric at those lambdas.
#'
#' @param object fitted \code{"cv.cpmnet"} object.
#' @param ... ignored.
#'
#' @return A list with selected lambdas and metrics, invisibly.
#'
#' @seealso \code{\link{cv.cpmnet}}, \code{\link{print.cv.cpmnet}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#' cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5, type = "brier",
#'                    nlambda = 10, parallel = FALSE)
#' summary(cvfit)
#' }
#'
#' @export
summary.cv.cpmnet <- function(object, ...) {
  if (!inherits(object, "cv.cpmnet")) {
    stop("object must inherit from class 'cv.cpmnet'", call. = FALSE)
  }
  cat("cv.cpmnet fit\n")
  cat(sprintf("  type             : %s\n", object$type))
  cat(sprintf("  nfolds_successful: %d\n", object$nfolds_successful))
  cat(sprintf("  nlambda          : %d\n", length(object$lambda)))
  if (object$type %in% c("mean", "median")) {
    cat(sprintf("  lambda.min (MAE) : %.6f\n", object$lambda.min.mae))
    cat(sprintf("  lambda.1se (MAE) : %.6f\n", object$lambda.1se.mae))
    cat(sprintf("  lambda.min (MSE) : %.6f\n", object$lambda.min.mse))
    cat(sprintf("  lambda.1se (MSE) : %.6f\n", object$lambda.1se.mse))
    if (!is.na(object$lambda.max.pR2)) {
      cat(sprintf("  lambda.max (pR2) : %.6f\n", object$lambda.max.pR2))
    }
    out <- list(
      lambda.min.mae = object$lambda.min.mae,
      lambda.1se.mae = object$lambda.1se.mae,
      lambda.min.mse = object$lambda.min.mse,
      lambda.1se.mse = object$lambda.1se.mse,
      lambda.max.pR2 = object$lambda.max.pR2
    )
  } else {
    cat(sprintf("  lambda.min       : %.6f\n", object$lambda.min))
    cat(sprintf("  lambda.1se       : %.6f\n", object$lambda.1se))
    idx_min <- which(object$lambda == object$lambda.min)
    idx_1se <- which(object$lambda == object$lambda.1se)
    if (length(idx_min) > 0L) {
      cat(sprintf("  cv metric @ min  : %.6f (se %.6f)\n",
                  object$cvm[idx_min], object$cvsd[idx_min]))
    }
    if (length(idx_1se) > 0L) {
      cat(sprintf("  cv metric @ 1se  : %.6f (se %.6f)\n",
                  object$cvm[idx_1se], object$cvsd[idx_1se]))
    }
    out <- list(
      lambda.min = object$lambda.min,
      lambda.1se = object$lambda.1se,
      cvm.min = if (length(idx_min) > 0L) object$cvm[idx_min] else NA_real_,
      cvm.1se = if (length(idx_1se) > 0L) object$cvm[idx_1se] else NA_real_
    )
  }
  invisible(out)
}
