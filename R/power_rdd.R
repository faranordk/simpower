#' Power and MDE for a sharp regression-discontinuity (RDD) design
#'
#' Simulation-based power analysis for the treatment jump at a cutoff in a sharp
#' RDD, estimated by local-linear regression with a triangular kernel. The null
#' is imposed by preserving the running-variable trend on each side while
#' permuting residuals (so there is no jump), and the RD estimator is refit on
#' each recombined sample. Because the local-linear estimator is linear in the
#' outcome, the analytic effect-size sweep is exact.
#'
#' Besides the overall power curve at the chosen bandwidth, a second result
#' traces how the MDE changes as the bandwidth changes.
#'
#' @inheritParams power_twfe
#' @param run Column name of the running (forcing) variable.
#' @param cutoff Threshold value of `run` at which treatment switches on
#'   (treatment = `run >= cutoff`).
#' @param bandwidth Half-width of the estimation window around the cutoff. When
#'   `NULL`, an MSE-optimal bandwidth from \pkg{rdrobust} is used if available,
#'   otherwise a rule-of-thumb value (and a message is shown).
#' @param bw_grid Optional vector of bandwidths for the MDE-versus-bandwidth
#'   graph. Defaults to a grid spanning roughly half to 1.5x the main bandwidth.
#' @param id Optional cluster id (SEs clustered by `id` when supplied).
#' @param bw_reps Number of simulations per bandwidth in the sweep (kept smaller
#'   than `reps` by default to bound run time).
#' @param kernel Weighting kernel for the local-linear fit: `"triangular"`
#'   (default) or `"uniform"` (equivalent to OLS within the bandwidth window).
#' @param n_units,n_treated Number of units and treated units (`run >= cutoff`)
#'   drawn from the *available* sample before the bandwidth window is applied.
#' @param use_rdrobust If `TRUE` and \pkg{rdrobust} is installed, use it for
#'   MSE-optimal bandwidth selection when `bandwidth` is `NULL`.
#'
#' @return An object of class `"simpower"` (design `"rdd"`). `extras$bandwidth`
#'   holds the MDE-versus-bandwidth table (its `n_in_bw` column counts
#'   *observations* inside each window) and `extras$h` the main bandwidth. In
#'   the sample-composition check, units and treated units are counted as
#'   clusters when `id` is supplied and as observations otherwise; the treated
#'   target refers to the full sample, the simulation average to what remains
#'   inside the bandwidth window.
#'
#' @examples
#' set.seed(1)
#' df <- sim_rdd(n = 1500, jump = 0)
#' pw <- power_rdd(df, y = "y", run = "x", cutoff = 0,
#'                 bandwidth = 0.5, reps = 200, bw_reps = 120)
#' pw
#' @export
power_rdd <- function(data, y, run, cutoff = 0, controls = NULL,
                      bandwidth = NULL, bw_grid = NULL, id = NULL,
                      n_units = NULL, n_treated = NULL,
                      reps = 1000, bw_reps = 300, emax = NULL, step = NULL,
                      kernel = c("triangular", "uniform"),
                      alpha = 0.05, alternative = c("greater", "two.sided", "less"),
                      crit_method = c("empirical", "normal"),
                      power_targets = c(0.8, 0.9), seed = 5000,
                      use_rdrobust = TRUE) {
  alternative <- match.arg(alternative)
  kernel <- match.arg(kernel)
  crit_method <- match.arg(crit_method)
  cl <- match.call()

  yv  <- .numcol(data, y, "y")
  rv  <- .numcol(data, run, "run")
  W   <- .controls_matrix(data, controls)
  idv <- if (!is.null(id)) .col(data, id, "id") else NULL

  keep0 <- stats::complete.cases(yv, rv, if (is.null(W)) NULL else W,
                                 if (is.null(idv)) NULL else idv)
  yv <- yv[keep0]; rv <- rv[keep0]
  if (!is.null(W)) W <- W[keep0, , drop = FALSE]
  run_c <- rv - cutoff
  treat <- as.numeric(run_c >= 0)
  n <- length(yv)
  uid <- if (!is.null(idv)) as.integer(factor(idv[keep0])) else seq_len(n)
  has_id <- !is.null(idv)

  ## bandwidth selection
  bwsel <- if (is.null(bandwidth)) .rdd_bw(yv, run_c, use_rdrobust) else
    list(h = bandwidth, source = "user-specified")
  h <- bwsel$h
  if (is.null(bandwidth)) message("simpower: RDD bandwidth = ", signif(h, 4),
                                  " (", bwsel$source, ")")

  treated_ids <- uid[treat == 1]
  real_units  <- length(unique(uid))
  real_treat  <- length(unique(treated_ids))
  sd_y <- stats::sd(yv[abs(run_c) <= h])                 # SD of Y within the window

  ## one simulation run at a given bandwidth
  rdd_run <- function(hh, R, sd) {
    one_rep <- function() {
      sm   <- .subsample_units(uid, treated_ids, n_units, n_treated)
      mask <- sm & (abs(run_c) <= hh)
      rc <- run_c[mask]; tr <- treat[mask]; y0 <- yv[mask]
      Wk <- if (is.null(W)) NULL else W[mask, , drop = FALSE]
      uk <- uid[mask]
      ## Diagnostics, all in consistent user-facing units:
      ##   units   -- clusters in the window when `id` is given, else obs
      ##   treated -- treated *units* on the same scale as `units` (and as the
      ##              `n_treated` target), not treated rows
      ##   nobs    -- raw observation count inside the window (for the sweep)
      dg <- function() c(
        units   = if (has_id) length(unique(uk)) else length(y0),
        treated = if (has_id) length(unique(uk[tr == 1])) else sum(tr == 1),
        nobs    = length(y0))
      if (length(unique(tr)) < 2L || length(y0) < 5L) {
        return(list(coef = c(effect = NA_real_), se = c(effect = NA_real_),
                    diag = dg()))
      }
      w  <- if (kernel == "uniform") rep(1, length(rc)) else pmax(0, 1 - abs(rc) / hh)
      Eexog <- cbind(rc, tr * rc, Wk)
      yk <- .recombine_resid(y0, Eexog)
      est <- .fit_rdd(yk, rc, tr, w, Wk, cluster = if (has_id) uk else NULL)
      list(coef = c(effect = unname(est["b"])),
           se   = c(effect = unname(est["se"])),
           diag = dg())
    }
    .run_sim(one_rep, R, sd)
  }

  ## main run at bandwidth h
  sim <- rdd_run(h, reps, seed)
  bs <- sim$B[, 1]; ses <- sim$S[, 1]
  .warn_dropped(bs, ses, reps)
  grid  <- auto_grid(ses, emax, step)
  power <- power_from_null(bs, ses, grid, alpha, alternative, crit_method)
  mde   <- mde_targets(grid, power, power_targets)

  ## MDE-versus-bandwidth sweep
  if (is.null(bw_grid)) {
    bw_grid <- sort(unique(round(h * seq(0.5, 1.5, by = 0.125), 6)))
  }
  th1 <- power_targets[1]
  sweep <- lapply(seq_along(bw_grid), function(i) {
    hh <- bw_grid[i]
    s  <- rdd_run(hh, bw_reps, seed + i)
    b  <- s$B[, 1]; se <- s$S[, 1]
    g  <- auto_grid(se, NULL, NULL)
    p  <- power_from_null(b, se, g, alpha, alternative, crit_method)
    m <- mde_from_curve(g, p, th1)
    data.frame(bandwidth = hh,
               n_in_bw   = mean(s$D[, "nobs"]),   # observations, as the name says
               se_mean   = mean(se, na.rm = TRUE),
               mde       = m,
               mde_std   = m / sd_y)
  })
  bandwidth_tbl <- do.call(rbind, sweep)
  attr(bandwidth_tbl, "target") <- th1

  ## Target and sim-average are now on the same scale: unique units (clusters
  ## when `id` is given, observations otherwise). Note the target counts treated
  ## units in the *full* sample while the sim average counts those left inside
  ## the bandwidth window -- the label makes that explicit.
  sample_check <- data.frame(
    metric  = c(if (has_id) "units, in bw (clusters)" else "obs in bandwidth",
                if (has_id) "treated units, in bw" else "treated obs, in bw"),
    target  = c(n_units %||% real_units, n_treated %||% real_treat),
    sim_avg = c(mean(sim$D[, "units"]), mean(sim$D[, "treated"]))
  )

  new_simpower(
    design = "rdd", call = cl,
    grid = grid, power = data.frame(effect = grid, power = power),
    mde = mde, bs = bs, ses = ses,
    alpha = alpha, alternative = alternative, reps = reps,
    request = list(n_units = n_units, n_treated = n_treated, controls = controls,
                   y = y, run = run, cutoff = cutoff, bandwidth = h),
    sample_check = sample_check,
    extras = list(h = h, bw_source = bwsel$source, bandwidth = bandwidth_tbl, sd_y = sd_y),
    data = data
  )
}

## Bandwidth selection: rdrobust MSE-optimal if available, else rule-of-thumb.
.rdd_bw <- function(y, run_c, use_rdrobust = TRUE) {
  if (use_rdrobust && requireNamespace("rdrobust", quietly = TRUE)) {
    bw <- tryCatch(rdrobust::rdbwselect(y, run_c)$bws[1, 1], error = function(e) NA_real_)
    if (is.finite(bw) && bw > 0) return(list(h = bw, source = "rdrobust (mserd)"))
  }
  h <- 0.5 * stats::IQR(run_c)
  if (!is.finite(h) || h <= 0) h <- stats::sd(run_c)
  list(h = h, source = "rule-of-thumb (0.5*IQR); install 'rdrobust' or set `bandwidth`")
}
