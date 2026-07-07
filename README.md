# cpmnet

Penalized Cumulative Probability Models with Elastic Net.

`cpmnet` fits cumulative probability models (CPMs) with an elastic-net
penalty along a full penalization (regularization) path, using coordinate descent
implemented in Fortran. Cross-validation selects `lambda` under one of
five conceptual criteria: mean or median regression loss, pinball
quantile loss, the Brier score over a percentile grid, or the held-out
log-likelihood. The Brier variant is the default and is the criterion
recommended in the accompanying methods paper.

## Installation

`cpmnet` depends on `rms` (for `orm.fit` and `probabilityFamilies`),
`MASS`, `PAmeasures`, and `parallel`. A Fortran compiler is required
at install time (macOS: `brew install gcc`; Debian/Ubuntu:
`apt install gfortran`; Windows: included in `Rtools`).

```r
# install.packages("devtools")
devtools::install_github("xuduovictorwen/cpmnet")
```

## Quick start

```r
library(cpmnet)

set.seed(1)
n <- 200
p <- 20
x <- matrix(rnorm(n * p), n, p)
beta_true <- c(1, -1, 0.5, rep(0, p - 3))
y <- x %*% beta_true + rnorm(n)

# fit the elastic net path (alpha_en = 0.5 is the default)
fit <- cpmnet(x, y, family = "probit", alpha_en = 0.5)

# one-page summary of the path and per-lambda diagnostics
summary(fit)

# coefficient-path plot (x-axis: log lambda / L1 norm / deviance explained)
plot(fit, xvar = "lambda")

# select lambda by 5-fold Brier-score CV (the default CV type)
cvfit <- cv.cpmnet(x, y, family = "probit", nfolds = 5,
                   type = "brier", parallel = FALSE)
plot(cvfit)

# extract slope coefficients at lambda.min
coef(cvfit, s = "lambda.min")

# include the ordinal intercepts as well
coef(cvfit, s = "lambda.min", include_alpha = TRUE)

# median prediction at lambda.min (type = "median" is the default)
predict(cvfit$cpmnet.fit, newx = x[1:5, ],
        type = "median", s = cvfit$lambda.min)

# mean prediction at lambda.min
predict(cvfit$cpmnet.fit, newx = x[1:5, ],
        type = "mean", s = cvfit$lambda.min)
```

## Arguments of interest

`cpmnet()`:

- `family`: one of `"logistic"`, `"probit"`, `"loglog"`, `"cloglog"`,
  `"cauchit"` (any name in `names(rms::probabilityFamilies)`).
- `alpha_en`: elastic-net mixing parameter in `[0, 1]` (`1` = lasso,
  `0` = ridge). The alternative name `alpha` is accepted for
  compatibility with glmnet call sites; `alpha_en` is the canonical
  name because `alpha` denotes the ordinal-intercept vector on the
  returned fit object.
- `weights`: observation weights. Zero-weight rows are dropped.
- `penalty_factor`: per-predictor penalty multipliers (length `p`).
- `adaptive = TRUE`: use inverse-|orm| weights for the adaptive lasso.
- `standardize = TRUE`: center and scale `x` before fitting (default).
- `verbose = TRUE`: warn when the path is truncated due to numerical
  issues; set to `FALSE` in simulation loops.

`cv.cpmnet()`:

- `type`: CV criterion, default `"brier"`. See below.
- `nfolds`, `foldid`: standard cross-validation controls.
- `parallel`, `ncores`: parallel folds via `parallel::mclapply`
  (non-Windows only; runs sequentially with a warning on Windows).
- `tau_levels`, `tau_weights`: quantile grid for `"pinball_*"`.
- `brier_probs`, `brier_weights`: percentile grid for `"brier"`.

## Cross-validation criteria

`cv.cpmnet(type = ...)` supports five `type` of conceptual
criteria emphasized in the methods paper:

- `"brier"` (default): squared Brier score averaged over response
  percentiles (deciles by default). Returns `lambda.min`, `lambda.1se`.
- `"mean"`, `"median"`: regression loss on predicted mean or median.
  Returns `lambda.min.mae`, `lambda.1se.mae`, `lambda.min.mse`,
  `lambda.1se.mse`, and (for `"mean"`) `lambda.max.pR2`,
  `lambda.1se.pR2`.
- `"pinball_abs"`, `"pinball_sq"`: pinball loss, absolute or squared,
  over `tau_levels`.
- `"loglik"`: held-out log-likelihood of the observed outcome
  categories, maximized (the criterion `ordinalNet` calls `cvLoglik`).
  `cvm` stores its negative, so `lambda.min` / `lambda.1se` keep their
  usual meaning. A held-out outcome outside the training support is
  assigned the boundary category.

## S3 methods

For a `cpmnet` fit: `coef`, `predict`, `summary`, `plot`, `print`.
For a `cv.cpmnet` fit: `coef` (with `s = "lambda.min"` or
`"lambda.1se"`), `summary`, `plot`, `print`.

## Testing

The package ships with a testthat suite of 131 tests (205 expectations)
covering input validation, all link families, all four CV types,
predict return shapes, unpenalized fits against `rms::orm.fit`,
intercept monotonicity across link functions, the adaptive lasso path,
and reproducibility. Run via `devtools::test()`.

## Citation

This package accompanies the manuscript *Penalized Cumulative
Probability Models* (forthcoming). Please cite the paper when using
`cpmnet` in research.

## Issues

Please report issues at
<https://github.com/xuduovictorwen/cpmnet/issues>.

## License

GPL (>= 2).
