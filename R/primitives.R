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
##
## `inject`: optional n-vector d, the direction a heterogeneous effect would be
## injected in (see the `het` option). Because the estimator is linear in y,
## fitting on y + e*d gives coefficient beta1(e) = beta1 + e*lam with
## lam = [(X'X)^-1 X'd]_1, and residuals u + e*uD with uD = d - X (X'X)^-1 X'd.
## The tested coefficient's sandwich variance is then EXACTLY quadratic in e:
## V11(e) = v0 + e*v1 + e^2*v2 (per-cluster scores are linear in e, the meat is
## their square). The returned `het = c(lam, v1, v2)` lets the power sweep
## evaluate the heterogeneous fit at every effect size without refitting.
.ols_vcov <- function(X, y, cluster = NULL, dfK = NULL, inject = NULL) {
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
  het <- if (!is.null(inject)) {
    .sandwich_inject(X, u, inject, XtXinv, cluster, fac)
  }
  list(beta = as.vector(beta), V = V, u = u, n = n, het = het)
}

## The heterogeneity pieces c(lam, v1, v2) for the FIRST regressor of a fitted
## sandwich (shared by OLS / 2SLS / weighted RDD): `Xbread` is the matrix whose
## cross-products form the meat scores (X, Xhat, or the weighted Xw), `u` the
## fitted residuals, `d` the injection direction ON THE SAME SCALE as u, `Binv`
## the bread inverse, `a` implied as Binv[1, ]. Exact by linearity.
.sandwich_inject <- function(Xbread, u, d, Binv, cluster, fac,
                             bD = NULL, Xresid = Xbread) {
  if (is.null(bD)) bD <- Binv %*% crossprod(Xbread, d)
  uD <- as.vector(d - Xresid %*% bD)       # residual response to the injection
  w  <- as.vector(Xbread %*% Binv[1L, ])   # row scores for coefficient 1
  if (is.null(cluster)) {
    s0 <- w * u; sD <- w * uD
  } else {
    cl <- as.integer(factor(cluster))
    s0 <- as.vector(rowsum(w * u, cl)); sD <- as.vector(rowsum(w * uD, cl))
  }
  c(lam = bD[1L],
    v1  = fac * 2 * sum(s0 * sD),
    v2  = fac * sum(sD * sD))
}

## Two-way fixed-effects estimate of the coefficient on the first regressor,
## with unit-clustered SE. `Xextra` holds optional control columns (already on
## the raw scale); FE are absorbed by within-demeaning. Returns c(b, se) for the
## tested regressor `xt`.
.fit_twfe <- function(y, xt, uid, ti, Xextra = NULL, inject = NULL) {
  G <- length(unique(uid)); Tn <- length(unique(ti))
  M <- cbind(xt, Xextra, y, inject)
  Md <- .demean_twoway(M, uid, ti)
  p  <- 1L + (if (is.null(Xextra)) 0L else ncol(as.matrix(Xextra)))
  Xd <- Md[, seq_len(p), drop = FALSE]
  yd <- Md[, p + 1L]
  dj <- if (is.null(inject)) NULL else Md[, p + 2L]
  dfK <- G + Tn + p                       # absorbed FE + slopes
  fit <- .ols_vcov(Xd, yd, cluster = uid, dfK = dfK, inject = dj)
  if (is.null(fit)) {
    out <- c(b = NA_real_, se = NA_real_)
    if (!is.null(inject)) out <- c(out, lam = NA_real_, v1 = NA_real_, v2 = NA_real_)
    return(out)
  }
  out <- c(b = fit$beta[1L], se = sqrt(max(fit$V[1L, 1L], 0)))
  if (!is.null(inject)) out <- c(out, fit$het)
  out
}

## Two-stage least squares. `D` is the (single) endogenous regressor tested,
## `W` the exogenous controls (WITHOUT intercept; intercept added here),
## `Zexcl` the excluded instrument(s). Coefficient on D is returned.
## Linear in y; robust (HC1) or cluster-robust SE.
.fit_tsls <- function(y, D, Zexcl, W = NULL, cluster = NULL, inject = NULL) {
  n <- length(y)
  bad <- function() {
    out <- c(b = NA_real_, se = NA_real_)
    if (!is.null(inject)) out <- c(out, lam = NA_real_, v1 = NA_real_, v2 = NA_real_)
    out
  }
  ones <- rep(1, n)
  X <- cbind(D, ones, W)                  # structural regressors (D first)
  Z <- cbind(Zexcl, ones, W)              # instruments + exogenous
  ZtZinv <- .safe_inv(crossprod(Z))
  if (is.null(ZtZinv)) return(bad())
  Xhat <- Z %*% (ZtZinv %*% crossprod(Z, X))
  A <- .safe_inv(crossprod(Xhat, X))
  if (is.null(A)) return(bad())
  beta <- A %*% crossprod(Xhat, y)
  u <- as.vector(y - X %*% beta)
  k <- ncol(X)
  if (n - k <= 0) return(bad())
  if (is.null(cluster)) {
    meat <- crossprod(Xhat * u); fac <- n / (n - k)
  } else {
    cl <- as.integer(factor(cluster)); G <- max(cl)
    if (G < 2L) return(bad())
    Sg <- rowsum(Xhat * u, cl); meat <- crossprod(Sg)
    fac <- (G / (G - 1)) * ((n - 1) / (n - k))
  }
  V <- A %*% meat %*% t(A) * fac
  out <- c(b = beta[1L], se = sqrt(max(V[1L, 1L], 0)))
  if (!is.null(inject)) {
    ## 2SLS is linear in y: beta(e) = beta + e * A Xhat'd, residuals against
    ## the STRUCTURAL X; meat scores use Xhat. Exact quadratic as elsewhere.
    bD <- A %*% crossprod(Xhat, inject)
    out <- c(out, .sandwich_inject(Xhat, u, inject, A, cluster, fac,
                                   bD = bD, Xresid = X))
  }
  out
}

## Weighted local-linear RDD estimate (sharp). `treat` = 1{run >= cutoff}.
## Model: y ~ treat + run_c + treat:run_c + W, weighted by kernel weights `w`.
## The coefficient on `treat` is the RD jump. Robust (HC1) or cluster SE.
.fit_rdd <- function(y, run_c, treat, w, W = NULL, cluster = NULL, inject = NULL) {
  n <- length(y)
  bad <- function() {
    out <- c(b = NA_real_, se = NA_real_)
    if (!is.null(inject)) out <- c(out, lam = NA_real_, v1 = NA_real_, v2 = NA_real_)
    out
  }
  ones <- rep(1, n)
  X <- cbind(treat, ones, run_c, treat * run_c, W)  # treat coef first
  sw <- sqrt(w)
  Xw <- X * sw; yw <- y * sw
  XtXinv <- .safe_inv(crossprod(Xw))
  if (is.null(XtXinv)) return(bad())
  beta <- XtXinv %*% crossprod(Xw, yw)
  u <- as.vector(yw - Xw %*% beta)
  k <- ncol(X)
  if (n - k <= 0) return(bad())
  if (is.null(cluster)) {
    meat <- crossprod(Xw * u); fac <- n / (n - k)
  } else {
    cl <- as.integer(factor(cluster)); G <- max(cl)
    if (G < 2L) return(bad())
    Sg <- rowsum(Xw * u, cl); meat <- crossprod(Sg)
    fac <- (G / (G - 1)) * ((n - 1) / (n - k))
  }
  V <- XtXinv %*% meat %*% XtXinv * fac
  out <- c(b = beta[1L], se = sqrt(max(V[1L, 1L], 0)))
  if (!is.null(inject)) {
    ## the injection lives on the raw outcome scale -> weight it like y
    out <- c(out, .sandwich_inject(Xw, u, inject * sw, XtXinv, cluster, fac))
  }
  out
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
  crit <- .null_crit(t0, alpha, alternative, crit_method)
  ## `grid` holds effect *magnitudes* (>= 0) in the hypothesised direction, so
  ## the injected true effect is -e for a one-sided "less" test.
  vapply(grid, function(e) {
    switch(alternative,
           greater   = mean((bs + e) / ses >  crit["hi"]),
           less      = mean((bs - e) / ses <  crit["lo"]),
           two.sided = mean(abs((bs + e) / ses) > crit["ab"]))
  }, numeric(1))
}

## Heterogeneous-effects power curve. `H` holds per-rep injection pieces
## (columns lam, v1, v2 -- see .ols_vcov). Injecting an average effect e with
## per-unit multipliers shifts the tested coefficient by e*lam and moves its
## sandwich variance to se^2 + e*v1 + e^2*v2, both EXACT by linearity -- so no
## refitting or bootstrapping across effect sizes is required. The critical
## value comes from the same simulated null t-statistics as the constant sweep
## (at e = 0 the two coincide), so the realized size is unchanged.
power_from_null_het <- function(bs, ses, H, grid, alpha = 0.05,
                                alternative = c("greater", "two.sided", "less"),
                                crit_method = c("empirical", "normal")) {
  alternative <- match.arg(alternative)
  crit_method <- match.arg(crit_method)
  ok <- is.finite(bs) & is.finite(ses) & ses > 0 &
        is.finite(H[, 1L]) & is.finite(H[, 2L]) & is.finite(H[, 3L])
  bs <- bs[ok]; ses <- ses[ok]
  lam <- H[ok, 1L]; v1 <- H[ok, 2L]; v2 <- H[ok, 3L]
  t0 <- bs / ses
  crit <- .null_crit(t0, alpha, alternative, crit_method)
  v0 <- ses^2
  vplus  <- function(e) sqrt(pmax(v0 + e * v1 + e^2 * v2, .Machine$double.eps))
  vminus <- function(e) sqrt(pmax(v0 - e * v1 + e^2 * v2, .Machine$double.eps))
  vapply(grid, function(e) {
    switch(alternative,
           greater   = mean((bs + e * lam) / vplus(e)  >  crit["hi"]),
           less      = mean((bs - e * lam) / vminus(e) <  crit["lo"]),
           two.sided = mean(abs((bs + e * lam) / vplus(e)) > crit["ab"]))
  }, numeric(1))
}

## Critical value(s) for the null t-statistics (shared by the constant and
## heterogeneous sweeps and the power-ceiling diagnostic).
.null_crit <- function(t0, alpha, alternative, crit_method) {
  if (crit_method == "normal") {
    q <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else stats::qnorm(1 - alpha)
    c(lo = -q, hi = q, ab = q)
  } else {
    c(lo = unname(stats::quantile(t0, alpha, na.rm = TRUE)),
      hi = unname(stats::quantile(t0, 1 - alpha, na.rm = TRUE)),
      ab = unname(stats::quantile(abs(t0), 1 - alpha, na.rm = TRUE)))
  }
}

## The power CEILING under heterogeneous effects. Because the injected
## heterogeneity inflates the estimated SE as the effect grows (the e^2 * v2
## term), the t-statistic has a finite asymptote lam / sqrt(v2) as e -> Inf,
## so power saturates below 1. If the ceiling sits below a requested power
## target, no average effect size -- however large -- can reach that target
## under the assumed heterogeneity; only more units can.
.het_power_ceiling <- function(bs, ses, H, alpha, alternative, crit_method) {
  ok <- is.finite(bs) & is.finite(ses) & ses > 0 &
        is.finite(H[, 1L]) & is.finite(H[, 3L])
  lam <- H[ok, 1L]; v2 <- H[ok, 3L]
  crit <- .null_crit(bs[ok] / ses[ok], alpha, alternative, crit_method)
  tinf <- ifelse(v2 > 0, lam / sqrt(v2), ifelse(lam > 0, Inf, ifelse(lam < 0, -Inf, 0)))
  switch(alternative,
         greater   = mean(tinf >  crit["hi"]),
         less      = mean(-tinf < crit["lo"]),
         two.sided = mean(abs(tinf) > crit["ab"]))
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
## `widen` stretches the default range (used with `het`, whose flatter curves
## reach the power targets at larger average effects).
auto_grid <- function(ses, emax = NULL, step = NULL, n_step = 200L,
                      widen = FALSE) {
  s <- stats::median(ses[is.finite(ses) & ses > 0])
  if (!is.finite(s) || s <= 0) s <- 1
  if (is.null(emax)) emax <- (if (widen) 8 else 5) * s
  if (is.null(step)) step <- emax / n_step
  seq(0, emax, by = step)
}

## Shared post-processing: the power curve and MDE, heterogeneity-aware when a
## `het` spec is active. With `het`, the returned curve/MDE reflect the
## heterogeneous injection, and `extras` carries the constant-effect
## counterpart (computed from the SAME simulated null draws -- common random
## numbers) so the cost of heterogeneity is directly readable.
.het_curves <- function(bs, ses, H, grid, alpha, alternative, crit_method,
                        power_targets, hspec, auto_extend = FALSE) {
  if (is.null(hspec)) {
    p <- power_from_null(bs, ses, grid, alpha, alternative, crit_method)
    return(list(power = p, mde = mde_targets(grid, p, power_targets),
                grid = grid, extras = NULL))
  }
  p    <- power_from_null_het(bs, ses, H, grid, alpha, alternative, crit_method)
  ceil <- .het_power_ceiling(bs, ses, H, alpha, alternative, crit_method)
  ## Heterogeneous curves rise more slowly, so an auto-chosen grid may stop
  ## short of the targets. When the user did not fix `emax` and the ceiling
  ## shows a target is attainable, stretch the grid (same number of points)
  ## until the reachable targets are on it.
  if (auto_extend) {
    reachable <- power_targets[power_targets <= ceil - 0.005]
    tries <- 0L
    while (tries < 6L && length(reachable) &&
           anyNA(mde_targets(grid, p, reachable, warn = FALSE))) {
      grid <- seq(0, max(grid) * 3, length.out = length(grid))
      p    <- power_from_null_het(bs, ses, H, grid, alpha, alternative, crit_method)
      tries <- tries + 1L
    }
  }
  m <- mde_targets(grid, p, power_targets, warn = FALSE)
  if (anyNA(m)) {
    miss   <- power_targets[is.na(m)]
    capped <- miss[miss > ceil - 0.005]
    if (length(capped)) {
      warning(sprintf(
        paste0("Under this effect heterogeneity, power plateaus at ~%.2f: the %s ",
               "power target(s) cannot be reached at ANY average effect size, ",
               "because the heterogeneity itself inflates the estimated SE as the ",
               "effect grows. A larger sample (not a larger effect) is needed."),
        ceil, paste0(capped * 100, "%", collapse = "/")), call. = FALSE)
    }
    if (length(setdiff(miss, capped))) {
      warning(sprintf(
        paste0("MDE for %s power not reached within emax = %.4g (max power %.2f, ",
               "ceiling %.2f). Increase `emax`."),
        paste0(setdiff(miss, capped) * 100, "%", collapse = "/"),
        max(grid), max(p, na.rm = TRUE), ceil), call. = FALSE)
    }
  }
  pc <- power_from_null(bs, ses, grid, alpha, alternative, crit_method)
  list(power = p, mde = m, grid = grid,
       extras = list(label = .het_label(hspec), spec = hspec$spec,
                     power_ceiling = ceil,
                     power_constant = data.frame(effect = grid, power = pc),
                     mde_constant = mde_targets(grid, pc, power_targets,
                                                warn = FALSE)))
}

.het_label <- function(hspec) {
  if (hspec$type == "multiset") {
    sprintf("multipliers {%s} (normalised to mean 1)",
            paste(signif(hspec$values, 3), collapse = ", "))
  } else "user-supplied multiplier function"
}
