#' cpmnet: Penalized Cumulative Probability Models with Elastic Net
#'
#' Fits cumulative probability models (CPMs) with elastic net
#' penalization (regularization) via coordinate descent, computed along
#' a full penalization (regularization) path. Cross-validation is provided for lambda
#' selection with four conceptual criteria: mean and median regression
#' loss, pinball quantile loss, and the Brier score; proper scoring
#' rule criteria are also available.
#' The core fitting loop is implemented in Fortran for computational
#' efficiency.
#'
#' @section Entry points:
#' \itemize{
#'   \item \code{\link{cpmnet}} fits a single model along a lambda path.
#'   \item \code{\link{cv.cpmnet}} performs k-fold cross-validation.
#'   \item \code{\link{predict.cpmnet}} produces point predictions.
#'   \item \code{\link{coef.cpmnet}}, \code{\link{summary.cpmnet}},
#'     \code{\link{plot.cpmnet}}, \code{\link{print.cpmnet}} provide
#'     S3 helpers for the fit object.
#' }
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
#' @useDynLib cpmnet, .registration = TRUE
#' @importFrom rms orm.fit probabilityFamilies
#' @importFrom MASS mvrnorm
#' @importFrom stats sd quantile
#' @importFrom utils head
#' @importFrom parallel detectCores mclapply
#' @importFrom grDevices rainbow
#' @importFrom graphics abline arrows axis grid legend matlines mtext points text plot
"_PACKAGE"
