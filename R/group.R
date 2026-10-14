# grouped penalization: a factor is penalized as one unit
# R centers and orthonormalizes each block, Fortran fits

# group labels to integer ids; 0 or "" means unpenalized, as in grpreg
#' @keywords internal
#' @noRd
normalize_group <- function(group, p) {
  if (is.null(group)) return(NULL)
  if (length(group) != p) {
    stop("length(group) = ", length(group), " must equal ncol(x) = ", p,
         call. = FALSE)
  }
  if (is.factor(group)) group <- as.character(group)
  unpen <- as.character(group) %in% c("0", "")
  lab <- as.character(group)
  gid <- integer(p)
  pen_lab <- unique(lab[!unpen])
  gid[!unpen] <- match(lab[!unpen], pen_lab)
  gid[unpen] <- 0L
  sizes <- if (length(pen_lab)) tabulate(gid[!unpen], nbins = length(pen_lab)) else integer(0)
  list(gid = gid, labels = pen_lab, sizes = sizes,
       n_unpenalized = sum(unpen),
       all_singleton = length(sizes) == 0L || max(sizes) == 1L)
}

# center every column and orthonormalize each block under the weighted inner
# product, so the penalty is on the block's fitted contribution and does not
# depend on which level was dropped. gamma is the largest eigenvalue of each
# block under the (1/n) inner product the core uses, 1 when the weights are equal
#' @keywords internal
#' @noRd
group_prep <- function(x, gid, wt) {
  n <- nrow(x)
  weighted <- !all(wt == wt[1])
  x_means <- if (weighted) colSums(x * wt) / sum(wt) else colMeans(x)
  xc <- sweep(x, 2, x_means)
  # unpenalized columns become blocks of their own
  unpen_cols <- which(gid == 0L)
  pen_ids <- setdiff(sort(unique(gid)), 0L)
  spec <- c(lapply(pen_ids, function(g) which(gid == g)),
            lapply(unpen_cols, function(j) j))
  pf0 <- c(rep(1, length(pen_ids)), rep(0, length(unpen_cols)))
  G <- length(spec)
  xt_blocks <- vector("list", G); bmats <- vector("list", G)
  cols_orig <- vector("list", G); dgs <- integer(G); gamma <- rep(1, G)
  for (g in seq_len(G)) {
    idx <- spec[[g]]
    xg <- xc[, idx, drop = FALSE]
    s <- if (weighted) crossprod(xg, xg * wt) / sum(wt) else crossprod(xg) / n
    ev <- eigen(s, symmetric = TRUE)
    keep <- ev$values > sqrt(.Machine$double.eps) * ev$values[1]
    q <- ev$vectors[, keep, drop = FALSE]; lam <- ev$values[keep]
    dh <- diag(1 / sqrt(lam), nrow = length(lam))
    xt_blocks[[g]] <- xg %*% q %*% dh
    bmats[[g]] <- q %*% dh
    cols_orig[[g]] <- idx; dgs[g] <- length(lam)
    if (weighted) {
      gamma[g] <- max(eigen(crossprod(xt_blocks[[g]]) / n, symmetric = TRUE,
                            only.values = TRUE)$values)
    }
  }
  xt <- do.call(cbind, xt_blocks)
  ends <- cumsum(dgs); starts <- ends - dgs + 1L
  list(xt = xt, ranges = Map(function(s, e) s:e, starts, ends),
       bmats = bmats, cols_orig = cols_orig, dgs = dgs, gamma = gamma,
       gstart = starts, x_means = x_means, G = G, pf0 = pf0)
}

# the grouped path, returned in the same shape as the elementwise one
#' @keywords internal
#' @noRd
cpmnet_group_path <- function(x, y, gid, family, alpha_en, lambda, nlambda,
                              lambda_min_ratio, lambda_max, penalty_factor,
                              weights, max_iter, tol, nulldev, verbose, mc) {
  yp <- cpm_prep_y(x, y, family, weights)
  x <- yp$x; wt <- yp$wt; y_int <- yp$y; k <- yp$k
  y_mapping <- yp$y_mapping; link <- yp$link; alpha_init <- yp$alpha_init
  n <- nrow(x); p <- ncol(x)

  prep <- group_prep(x, gid, wt)
  G <- prep$G
  # user factors cover the penalized groups; unpenalized blocks come after them
  pf <- if (is.null(penalty_factor)) prep$pf0
        else c(penalty_factor, rep(0, sum(prep$pf0 == 0)))
  if (length(pf) != G) {
    stop("penalty_factor must have length = number of groups (",
         sum(prep$pf0 > 0), ")", call. = FALSE)
  }

  ptilde <- ncol(prep$xt)
  ni <- as.integer(n); ki <- as.integer(k); pti <- as.integer(ptilde)
  gi <- as.integer(G); linki <- as.integer(link)
  yf <- as.integer(y_int); wtf <- as.double(wt)
  xtf <- matrix(as.double(prep$xt), n, ptilde)
  gstart <- as.integer(prep$gstart); gsize <- as.integer(prep$dgs)
  pff <- as.double(pf); gamf <- as.double(prep$gamma); aif <- as.double(alpha_init)

  first_null <- FALSE
  if (is.null(lambda)) {
    if (!is.null(lambda_max)) {
      lmax <- lambda_max
    } else {
      # at a computed lambda_max every penalized group is zero unless alpha_en is floored
      first_null <- alpha_en >= .CPMNET_LAMBDA_MAX_ALPHA_FLOOR
      # an unpenalized group must not set lambda_max; a huge factor removes it
      pf_lmax <- pff
      pf_lmax[pf_lmax <= 0] <- 1e12
      lm <- .Fortran("cpmnet2_lambda_max_group", ni, ki, pti, gi, xtf, yf, wtf,
                     linki, gstart, gsize, pf_lmax, aif,
                     lmax = double(1), salloc = integer(1), PACKAGE = "cpmnet")
      if (lm$salloc != 0) stop("Fortran allocation failed in lambda_max", call. = FALSE)
      if (!is.finite(lm$lmax) || lm$lmax <= 0) {
        stop("lambda_max is not finite; every group may be unpenalized", call. = FALSE)
      }
      lmax <- lm$lmax / max(alpha_en, .CPMNET_LAMBDA_MAX_ALPHA_FLOOR)
      # lambda_max is taken with the unpenalized groups at zero rather than at
      # their optimum, so raise it until every penalized group is out
      if (any(pf <= 0)) {
        pen_rows <- unlist(prep$ranges[pf > 0])
        for (tries in 1:20) {
          probe <- .Fortran("cpmnet2_path_group", ni, ki, pti, gi, xtf, yf, wtf,
                            linki, gstart, gsize, pff, gamf, 1L, as.double(lmax),
                            as.double(alpha_en), aif, double(ptilde),
                            as.integer(max_iter), as.double(tol),
                            alpha_mat = matrix(0.0, k, 1L),
                            beta_mat = matrix(0.0, ptilde, 1L),
                            logL = double(1), converged = logical(1),
                            niter = integer(1), salloc = integer(1),
                            shuffle_seed = .cpmnet_shuffle_seed(),
                            PACKAGE = "cpmnet")
          if (probe$salloc != 0) break
          # a rounding leftover near 1e-15 counts as zero
          if (all(abs(probe$beta_mat[pen_rows, 1]) < 1e-10)) break
          lmax <- lmax * 1.5
        }
      }
    }
    if (is.null(lambda_min_ratio))
      lambda_min_ratio <- default_lambda_min_ratio(n, p, family)
    lambda <- create_lambda_path(lmax, lmax * lambda_min_ratio, nlambda)
  } else {
    lambda <- sort(lambda, decreasing = TRUE)
  }
  nlambda <- length(lambda)
  # the first fit takes a far larger penalty so the penalized groups are exactly zero
  # rather than near 1e-15; unpenalized groups have no penalty at any lambda
  lambda_fit <- lambda
  if (first_null) lambda_fit[1] <- 1e10 * lambda[1]

  fit <- .Fortran("cpmnet2_path_group", ni, ki, pti, gi, xtf, yf, wtf, linki,
                  gstart, gsize, pff, gamf, as.integer(nlambda), as.double(lambda_fit),
                  as.double(alpha_en), aif, double(ptilde),
                  as.integer(max_iter), as.double(tol),
                  alpha_mat = matrix(0.0, k, nlambda),
                  beta_mat = matrix(0.0, ptilde, nlambda),
                  logL = double(nlambda), converged = logical(nlambda),
                  niter = integer(nlambda), salloc = integer(1),
                  shuffle_seed = .cpmnet_shuffle_seed(),
                  PACKAGE = "cpmnet")
  st <- .path_status(fit, lambda, max_iter, verbose, "cpmnet2_path_group")
  fit <- st$res
  lambda <- st$lambda
  nlambda <- length(lambda)

  # back to the original columns, and the count of nonzero groups
  beta_orig <- matrix(0, p, nlambda); ngrp <- integer(nlambda)
  for (il in seq_len(nlambda)) {
    ng <- 0L
    for (g in seq_len(G)) {
      btg <- fit$beta_mat[prep$ranges[[g]], il]
      beta_orig[prep$cols_orig[[g]], il] <- as.numeric(prep$bmats[[g]] %*% btg)
      if (sqrt(sum(btg^2)) > 0) ng <- ng + 1L
    }
    ngrp[il] <- ng
  }
  adj <- as.numeric(crossprod(prep$x_means, beta_orig))
  alpha_orig <- fit$alpha_mat - matrix(adj, k, nlambda, byrow = TRUE)

  if (is.null(nulldev)) {
    nulldev <- rms::orm.fit(x = NULL, y = y_mapping[y_int + 1],
                            family = family, weights = wt)$deviance
  }

  dp <- list(x = x, y = y_int, wt = wt, n = n, p = p, k = k, link = link,
             alpha_init = alpha_init, beta_init = rep(0, p),
             x_means = prep$x_means, x_sds = rep(1, p), y_mapping = y_mapping)

  structure(list(
    alpha = alpha_orig, beta = beta_orig,
    alpha_scaled = fit$alpha_mat, beta_scaled = fit$beta_mat,
    df = colSums(beta_orig != 0), ngroups = ngrp,
    dim = c(p, nlambda), lambda = lambda,
    dev.ratio = pmax(0, 1 - fit$logL / nulldev), nulldev = nulldev,
    logL = fit$logL, converged = fit$converged, niter = fit$niter,
    nobs = n, data_prep = dp, alpha_en = alpha_en, standardize = TRUE,
    group = gid, group_labels = attr(gid, "labels"),
    beta_tilde = fit$beta_mat, call = mc
  ), class = "cpmnet")
}
