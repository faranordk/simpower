## ---------------------------------------------------------------------------
## Low-level statistical primitives (internal, not exported)
##
## All estimators here return the tested coefficient(s) together with a
## standard error, and every one of them is *linear in the outcome y*. That
## linearity is what makes the analytic effect-size sweep exact (see
## power_from_null()).
## ---------------------------------------------------------------------------

## Group-demean the columns of a matrix by an integer group index in 1..G.
.gdemean <- function(M, g) {
  gs   <- rowsum(M, g)              # G x k group sums (rows ordered by group id)
  cnt  <- tabulate(g)               # counts per group id
  means <- gs / cnt
  M - means[g, , drop = FALSE]
}

## Two-way within transformation via alternating projections.
## For a balanced panel one pass is exact; for unbalanced panels the two
## demeanings do not commute, so we iterate to a tolerance --- this reproduces
## what reghdfe / feols get by absorbing both fixed effects.
.demean_twoway <- function(M, uid, ti, tol = 1e-10, maxit = 200L) {
  M <- as.matrix(M)
  ## For a balanced panel the two demeanings commute, so a single pass of each is
  ## exact -- skip the iteration entirely. Balance <=> every unit has the same
  ## number of periods and every period the same number of units.
  Tn <- length(unique(ti)); G <- length(unique(uid))
  if (all(tabulate(uid) == Tn) && all(tabulate(ti) == G)) {
    return(.gdemean(.gdemean(M, uid), ti))
  }
  for (it in seq_len(maxit)) {
    M0 <- M
    M  <- .gdemean(M, uid)
    M  <- .gdemean(M, ti)
    if (max(abs(M - M0)) < tol) break
  }
  M
}

## Safe inverse of a small cross-product matrix.
.safe_inv <- function(A) {
  out <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (is.null(out)) out <- tryCatch(solve(A), error = function(e) NULL)
  out
}

## OLS with (optional) cluster-robust or HC1 standard errors.
## Returns the coefficient vector, its variance matrix, and residuals.
## dfK is the number of estimated parameters used in the small-sample factor
## (for absorbed fixed effects this includes the FE count).
.ols_vcov <- function(X, y, cluster = NULL, dfK = NULL) {
  X <- as.matrix(X)
  n <- nrow(X); k <- ncol(X)
  XtXinv <- .safe_inv(crossprod(X))
  if (is.null(XtXinv)) return(NULL)
  beta <- XtXinv %*% crossprod(X, y)
  u    <- as.vector(y - X %*% beta)
  if (is.null(dfK)) dfK <- k
  ## Too few residual degrees of freedom -> the small-sample factor would be
  ## non-positive and the variance meaningless. Signal an invalid fit (NULL) so
  ## the caller drops this rep instead of silently reporting SE = 0.
  if (n - dfK <= 0) return(NULL)
  if (is.null(cluster)) {                                   # HC1
    meat <- crossprod(X * u)
    fac  <- n / (n - dfK)
  } else {                                                  # CR1 (Stata-style)
    cl <- as.integer(factor(cluster)); G <- max(cl)
    if (G < 2L) return(NULL)                                # need >= 2 clusters
    Sg <- rowsum(X * u, cl)
    meat <- crossprod(Sg)
    fac  <- (G / (G - 1)) * ((n - 1) / (n - dfK))
  }
  V <- XtXinv %*% meat %*% XtXinv * fac
  list(beta = as.vector(beta), V = V, u = u, n = n)
}

## Two-way fixed-effects estimate of the coefficient on the first regressor,
## with unit-clustered SE. `Xextra` holds optional control columns (already on
## the raw scale); FE are absorbed by within-demeaning. Returns c(b, se) for the
## tested regressor `xt`.
.fit_twfe <- function(y, xt, uid, ti, Xextra = NULL) {
  G <- length(unique(uid)); Tn <- length(unique(ti))
  M <- cbind(xt, Xextra, y)
  Md <- .demean_twoway(M, uid, ti)
  p  <- 1L + (if (is.null(Xextra)) 0L else ncol(as.matrix(Xextra)))
  Xd <- Md[, seq_len(p), drop = FALSE]
  yd <- Md[, p + 1L]
  dfK <- G + Tn + p                       # absorbed FE + slopes
  fit <- .ols_vcov(Xd, yd, cluster = uid, dfK = dfK)
  if (is.null(fit)) return(c(b = NA_real_, se = NA_real_))
  c(b = fit$beta[1L], se = sqrt(max(fit$V[1L, 1L], 0)))
}

## Two-stage least squares. `D` is the (single) endogenous regressor tested,
## `W` the exogenous controls (WITHOUT intercept; intercept added here),
## `Zexcl` the excluded instrument(s). Coefficient on D is returned.
## Linear in y; robust (HC1) or cluster-robust SE.
.fit_tsls <- function(y, D, Zexcl, W = NULL, cluster = NULL) {
  n <- length(y)
  ones <- rep(1, n)
  X <- cbind(D, ones, W)                  # structural regressors (D first)
  Z <- cbind(Zexcl, ones, W)              # instruments + exogenous
  ZtZinv <- .safe_inv(crossprod(Z))
  if (is.null(ZtZinv)) return(c(b = NA_real_, se = NA_real_))
  Xhat <- Z %*% (ZtZinv %*% crossprod(Z, X))
  A <- .safe_inv(crossprod(Xhat, X))
  if (is.null(A)) return(c(b = NA_real_, se = NA_real_))
  beta <- A %*% crossprod(Xhat, y)
  u <- as.vector(y - X %*% beta)
  k <- ncol(X)
  if (n - k <= 0) return(c(b = NA_real_, se = NA_real_))
  if (is.null(cluster)) {
    meat <- crossprod(Xhat * u); fac <- n / (n - k)
  } else {
    cl <- as.integer(factor(cluster)); G <- max(cl)
    if (G < 2L) return(c(b = NA_real_, se = NA_real_))
    Sg <- rowsum(Xhat * u, cl); meat <- crossprod(Sg)
    fac <- (G / (G - 1)) * ((n - 1) / (n - k))
  }
  V <- A %*% meat %*% t(A) * fac
  c(b = beta[1L], se = sqrt(max(V[1L, 1L], 0)))
}

## Weighted local-linear RDD estimate (sharp). `treat` = 1{run >= cutoff}.
## Model: y ~ treat + run_c + treat:run_c + W, weighted by kernel weights `w`.
## The coefficient on `treat` is the RD jump. Robust (HC1) or cluster SE.
.fit_rdd <- function(y, run_c, treat, w, W = NULL, cluster = NULL) {
  n <- length(y)
  ones <- rep(1, n)
  X <- cbind(treat, ones, run_c, treat * run_c, W)  # treat coef first
  sw <- sqrt(w)
  Xw <- X * sw; yw <- y * sw
  XtXinv <- .safe_inv(crossprod(Xw))
  if (is.null(XtXinv)) return(c(b = NA_real_, se = NA_real_))
  beta <- XtXinv %*% crossprod(Xw, yw)
  u <- as.vector(yw - Xw %*% beta)
  k <- ncol(X)
  if (n - k <= 0) return(c(b = NA_real_, se = NA_real_))
  if (is.null(cluster)) {
    meat <- crossprod(Xw * u); fac <- n / (n - k)
  } else {
    cl <- as.integer(factor(cluster)); G <- max(cl)
    if (G < 2L) return(c(b = NA_real_, se = NA_real_))
    Sg <- rowsum(Xw * u, cl); meat <- crossprod(Sg)
    fac <- (G / (G - 1)) * ((n - 1) / (n - k))
  }
  V <- XtXinv %*% meat %*% XtXinv * fac
  c(b = beta[1L], se = sqrt(max(V[1L, 1L], 0)))
}

## ---------------------------------------------------------------------------
## Power curve and MDE from the stored null draws.
##
## Adding a constant true effect e shifts every stored coefficient by exactly e
## and leaves the SE unchanged, so power(e) = share of injected t-stats that
## clear the critical value. This is the analytic sweep (Doucette 2025, sec. 3).
## ---------------------------------------------------------------------------
power_from_null <- function(bs, ses, grid, alpha = 0.05,
                            alternative = c("greater", "two.sided", "less"),
                            crit_method = c("empirical", "normal")) {
  alternative <- match.arg(alternative)
  crit_method <- match.arg(crit_method)
  ok <- is.finite(bs) & is.finite(ses) & ses > 0
  bs <- bs[ok]; ses <- ses[ok]
  ## Critical value for the t-statistic under the imposed null.
  ##
  ## "normal" uses the closed-form Gaussian quantile. "empirical" (default)
  ## reads the critical value off the simulated null distribution of the
  ## t-statistics themselves, so the test has *exactly* the requested size even
  ## when clustered / small-sample SEs make the null heavier-tailed than normal.
  ## This is the more faithful "use the design's own behaviour" choice and it is
  ## free -- the null draws are already in hand.
  t0 <- bs / ses
  crit <- if (crit_method == "normal") {
    q <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else stats::qnorm(1 - alpha)
    c(lo = -q, hi = q, ab = q)
  } else {
    c(lo = unname(stats::quantile(t0, alpha, na.rm = TRUE)),
      hi = unname(stats::quantile(t0, 1 - alpha, na.rm = TRUE)),
      ab = unname(stats::quantile(abs(t0), 1 - alpha, na.rm = TRUE)))
  }
  ## `grid` holds effect *magnitudes* (>= 0) in the hypothesised direction, so
  ## the injected true effect is -e for a one-sided "less" test.
  vapply(grid, function(e) {
    switch(alternative,
           greater   = mean((bs + e) / ses >  crit["hi"]),
           less      = mean((bs - e) / ses <  crit["lo"]),
           two.sided = mean(abs((bs + e) / ses) > crit["ab"]))
  }, numeric(1))
}

## First grid point where the power curve reaches `th`.
mde_from_curve <- function(grid, power, th) {
  i <- which(power >= th)
  if (length(i)) grid[i[1L]] else NA_real_
}

## MDE at each requested power target, with a warning when a target is not
## reached within the effect grid (i.e. `emax` is too small and MDE is NA).
mde_targets <- function(grid, power, targets, warn = TRUE) {
  m <- stats::setNames(
    vapply(targets, function(th) mde_from_curve(grid, power, th), numeric(1)),
    paste0(targets * 100, "%"))
  if (warn && anyNA(m)) {
    miss <- targets[is.na(m)]
    warning(sprintf(
      paste0("MDE for %s power not reached within emax = %.4g ",
             "(max power attained = %.2f). Increase `emax` (or leave it NULL to auto-scale)."),
      paste0(miss * 100, "%", collapse = "/"), max(grid), max(power, na.rm = TRUE)),
      call. = FALSE)
  }
  m
}

## Auto grid on the effect scale when the user does not supply emax/step.
auto_grid <- function(ses, emax = NULL, step = NULL, n_step = 200L) {
  s <- stats::median(ses[is.finite(ses) & ses > 0])
  if (!is.finite(s) || s <= 0) s <- 1
  if (is.null(emax)) emax <- 5 * s
  if (is.null(step)) step <- emax / n_step
  seq(0, emax, by = step)
}
