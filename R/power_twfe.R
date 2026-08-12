#' Power and MDE for a two-way fixed-effects (TWFE) design
#'
#' Simulation-based power analysis for the coefficient \eqn{\beta} on a
#' treatment/independent variable in a two-way (unit + time) fixed-effects
#' model with unit-clustered standard errors:
#' \deqn{Y_{it} = \beta X_{it} + \gamma' Z_{it} + \alpha_i + \delta_t + \varepsilon_{it}.}
#'
#' The outcome is recombined across units to impose the empirical null (see
#' [simpower-package]); the design's own within-estimator and cluster-robust SE
#' are computed on each recombined sample; and a growing true effect is injected
#' analytically to trace the power curve.
#'
#' @param data A data frame in long (panel) format.
#' @param y Column name (string) of the outcome variable.
#' @param x Column name of the treatment / independent variable whose
#'   coefficient is tested.
#' @param id Column name of the unit id (panel variable).
#' @param time Column name of the time variable.
#' @param controls Optional character vector of control-variable column names.
#'   Factors/characters are expanded to dummies. Controls are partialled out in
#'   every simulated fit.
#' @param n_units Optional. Total number of units to draw for each simulated
#'   sample. When `NULL` the full sample is used.
#' @param n_treated Optional. Number of treated units (units with within-unit
#'   variation in `x`) to draw per simulated sample. When `NULL` (and `n_units`
#'   is given) the observed treated share is preserved.
#' @param reps Number of simulations (curve smoothness; does not affect bias).
#' @param emax,step Optional largest effect size and grid step on the power
#'   curve. When `NULL` they are chosen automatically from the null SE scale.
#' @param alpha Significance level.
#' @param alternative Test direction: `"greater"` (one-sided, the default and
#'   the paper's convention), `"less"`, or `"two.sided"`.
#' @param crit_method How the critical value for the test is obtained.
#'   `"empirical"` (default) reads it off the simulated null distribution of the
#'   t-statistics, giving the test exactly the requested size even under
#'   clustered / small-sample standard errors. `"normal"` uses the closed-form
#'   Gaussian quantile (the previous behaviour).
#' @param power_targets Power levels at which to report the MDE.
#' @param seed Random seed (for reproducibility).
#' @param use_fixest If `TRUE` and the \pkg{fixest} package is installed, use
#'   `fixest::feols()` for the real-data point estimate reported in
#'   `extras$estimate`; otherwise the internal within-estimator (which matches
#'   `feols`) is used. The simulation itself always uses the internal estimator.
#'
#' @return An object of class `"simpower"`. See [summary.simpower()] and
#'   [plot.simpower()].
#'
#' @examples
#' set.seed(1)
#' df <- sim_panel(n_units = 60, n_time = 8)
#' pw <- power_twfe(df, y = "y", x = "x", id = "id", time = "t", reps = 200)
#' pw
#' @export
power_twfe <- function(data, y, x, id, time, controls = NULL,
                       n_units = NULL, n_treated = NULL,
                       reps = 1000, emax = NULL, step = NULL,
                       alpha = 0.05, alternative = c("greater", "two.sided", "less"),
                       crit_method = c("empirical", "normal"),
                       power_targets = c(0.8, 0.9), seed = 5000,
                       use_fixest = TRUE) {
  alternative <- match.arg(alternative)
  crit_method <- match.arg(crit_method)
  cl <- match.call()

  yv <- .numcol(data, y, "y")
  xv <- .numcol(data, x, "x")
  idv <- .col(data, id, "id")
  tv  <- .col(data, time, "time")
  W   <- .controls_matrix(data, controls)

  keep0 <- stats::complete.cases(yv, xv, idv, tv, if (is.null(W)) NULL else W)
  yv <- yv[keep0]; xv <- xv[keep0]
  uid_full <- as.integer(factor(idv[keep0]))
  ti_full  <- as.integer(factor(tv[keep0]))
  if (!is.null(W)) W <- W[keep0, , drop = FALSE]

  ## treated units = those with within-unit variation in x
  rng0 <- tapply(xv, uid_full, function(v) diff(range(v)))
  treated_ids <- as.integer(names(rng0))[rng0 > 0]
  real_units  <- length(unique(uid_full))
  real_treat  <- length(treated_ids)
  sd_y <- stats::sd(yv - stats::ave(yv, uid_full))       # within-unit SD of Y

  one_rep <- function() {
    sm <- .subsample_units(uid_full, treated_ids, n_units, n_treated)
    u  <- as.integer(factor(uid_full[sm]))
    t  <- as.integer(factor(ti_full[sm]))
    xs <- xv[sm]; ys <- yv[sm]
    Ws <- if (is.null(W)) NULL else W[sm, , drop = FALSE]

    rec  <- .recombine_panel(ys, u, t)
    keep <- rec$keep
    yk <- rec$y[keep]
    uk <- as.integer(factor(u[keep])); tk <- as.integer(factor(t[keep]))
    xk <- xs[keep]
    Wk <- if (is.null(Ws)) NULL else Ws[keep, , drop = FALSE]

    est <- .fit_twfe(yk, xk, uk, tk, Wk)
    rng <- tapply(xk, uk, function(v) diff(range(v)))
    list(coef = c(effect = unname(est["b"])),
         se   = c(effect = unname(est["se"])),
         diag = c(units = length(unique(uk)), treated = sum(rng > 0)))
  }

  sim <- .run_sim(one_rep, reps, seed)
  bs <- sim$B[, 1]; ses <- sim$S[, 1]
  .warn_dropped(bs, ses, reps)
  grid  <- auto_grid(ses, emax, step)
  power <- power_from_null(bs, ses, grid, alpha, alternative, crit_method)
  mde   <- mde_targets(grid, power, power_targets)

  ## real-data point estimate + clustered SE (diagnostic only)
  ref <- .twfe_reference(yv, xv, uid_full, ti_full, W, use_fixest)

  req_units   <- n_units   %||% real_units
  req_treated <- n_treated %||% real_treat
  sample_check <- data.frame(
    metric  = c("units (clusters)", "treated units"),
    target  = c(req_units, req_treated),
    sim_avg = c(mean(sim$D[, "units"]), mean(sim$D[, "treated"]))
  )

  new_simpower(
    design = "twfe", call = cl,
    grid = grid, power = data.frame(effect = grid, power = power),
    mde = mde, bs = bs, ses = ses,
    alpha = alpha, alternative = alternative, reps = reps,
    request = list(n_units = n_units, n_treated = n_treated, controls = controls,
                   y = y, x = x, id = id, time = time),
    sample_check = sample_check,
    extras = list(use_fixest = use_fixest, sd_y = sd_y,
                  estimate = ref$estimate, estimator = ref$source),
    data = data
  )
}

## Real-data TWFE point estimate on the tested regressor, with unit-clustered SE.
## Uses fixest::feols() when `use_fixest` is TRUE and the package is available;
## otherwise (and identically, up to floating point) the internal within
## estimator. Wrapped so a failure only costs the diagnostic, not the fit.
.twfe_reference <- function(y, x, uid, ti, W = NULL, use_fixest = TRUE) {
  if (use_fixest && requireNamespace("fixest", quietly = TRUE)) {
    est <- tryCatch({
      df <- data.frame(y = y, x = x, .id = factor(uid), .t = factor(ti))
      rhs <- "x"
      if (!is.null(W)) { for (j in seq_len(ncol(W))) df[[paste0("w", j)]] <- W[, j]
                         rhs <- paste(c("x", paste0("w", seq_len(ncol(W)))), collapse = " + ") }
      fml <- stats::as.formula(sprintf("y ~ %s | .id + .t", rhs))
      m  <- fixest::feols(fml, data = df, cluster = ~ .id)
      co <- summary(m)$coeftable
      c(b = unname(co["x", 1L]), se = unname(co["x", 2L]))
    }, error = function(e) NULL)
    if (!is.null(est)) return(list(estimate = est, source = "fixest::feols"))
  }
  est <- .fit_twfe(y, x, uid, ti, W)
  list(estimate = c(b = unname(est["b"]), se = unname(est["se"])),
       source = "internal within-estimator")
}
