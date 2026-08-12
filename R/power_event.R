#' Power and MDE for an event-study design
#'
#' Simulation-based power analysis for a dynamic (event-study) treatment effect
#' in a two-way fixed-effects panel. Two objects are produced:
#'
#' * an **overall** power curve for a single shaped-exposure coefficient. The
#'   assumed effect path across post-treatment horizons is set by `shape`
#'   (`"constant"`, `"increasing"`, or `"decreasing"`), normalised so the peak
#'   horizon has weight 1; the overall MDE is therefore expressed in units of
#'   the **peak** effect.
#' * a **per-horizon** minimum detectable effect: for each post-treatment
#'   horizon \eqn{h}, the smallest effect the design could detect at that
#'   horizon on its own, read off the dynamic event-study regression. This is
#'   reported both as the raw effect (`mde`) and, taking the assumed `shape`
#'   into account, in units of the peak effect (`mde_peak = mde / weight`); the
#'   per-horizon plot uses the latter.
#'
#' Treatment timing is taken from an absorbing 0/1 indicator `treat` (once on it
#' stays on) or, alternatively, from a `cohort` column giving each unit's
#' adoption time. Relative event time is `time - adoption`; horizon `-1` is the
#' omitted reference period.
#'
#' @inheritParams power_twfe
#' @param treat Column name of an absorbing 0/1 treatment indicator. Optional if
#'   `cohort` is supplied.
#' @param cohort Optional column name giving each unit's adoption time (on the
#'   same scale as `time`); never-treated units take `NA`/`Inf`.
#' @param shape Assumed shape of the true effect path over post-treatment
#'   horizons: `"constant"`, `"increasing"`, or `"decreasing"`.
#' @param horizons Maximum post-treatment horizon to include -- i.e. how many
#'   periods after treatment to model (default: as many as the data support,
#'   capped at 20).
#' @param leads Number of pre-treatment leads to include as placebo terms
#'   (default: up to 5). Horizon `-1` is always the omitted baseline.
#' @param n_treated Number of eventually-treated units to draw per simulated
#'   sample (see [power_twfe()]).
#'
#' @return An object of class `"simpower"` (design `"event"`). Its `extras$horizon`
#'   element holds the per-horizon MDE table.
#'
#' @examples
#' set.seed(1)
#' df <- sim_event(n_units = 80, n_time = 10, adopt = 6)
#' pw <- power_event(df, y = "y", id = "id", time = "t", treat = "d",
#'                   shape = "increasing", reps = 150)
#' pw
#' @export
power_event <- function(data, y, id, time, treat = NULL, cohort = NULL,
                        controls = NULL,
                        shape = c("constant", "increasing", "decreasing"),
                        horizons = NULL, leads = NULL,
                        n_units = NULL, n_treated = NULL,
                        reps = 1000, emax = NULL, step = NULL,
                        alpha = 0.05, alternative = c("greater", "two.sided", "less"),
                        crit_method = c("empirical", "normal"),
                        power_targets = c(0.8, 0.9), seed = 5000) {
  shape <- match.arg(shape)
  alternative <- match.arg(alternative)
  crit_method <- match.arg(crit_method)
  cl <- match.call()

  yv  <- .numcol(data, y, "y")
  idv <- .col(data, id, "id")
  tv  <- .col(data, time, "time")
  W   <- .controls_matrix(data, controls)
  if (is.null(treat) && is.null(cohort)) {
    stop("Supply either `treat` (absorbing 0/1) or `cohort` (adoption time).",
         call. = FALSE)
  }

  keep0 <- stats::complete.cases(yv, idv, tv, if (is.null(W)) NULL else W)
  yv <- yv[keep0]
  uid <- as.integer(factor(idv[keep0]))
  ## integer time index preserving order of the original time values
  tlev <- sort(unique(tv[keep0]))
  ti  <- match(tv[keep0], tlev)
  if (!is.null(W)) W <- W[keep0, , drop = FALSE]
  N <- max(uid)

  ## adoption time (as a time index) per unit
  if (!is.null(cohort)) {
    cohv <- .col(data, cohort, "cohort")[keep0]
    g_idx <- match(cohv, tlev)                    # NA if never-treated / out of range
    gvec  <- tapply(g_idx, uid, function(v) v[1L])
  } else {
    trv <- .numcol(data, treat, "treat")[keep0]
    gvec <- tapply(seq_along(uid), uid, function(ix) {
      on <- ix[trv[ix] == 1]
      if (length(on)) min(ti[on]) else NA_integer_
    })
  }
  g_by_row <- gvec[uid]
  rel <- ti - g_by_row                            # relative event time; NA if never-treated

  ## horizon window
  post_avail <- sort(unique(rel[is.finite(rel) & rel >= 0]))
  if (!length(post_avail)) stop("No post-treatment periods found.", call. = FALSE)
  H <- if (is.null(horizons)) min(max(post_avail), 20L) else min(horizons, max(post_avail))
  Kpost <- intersect(0:H, post_avail)

  pre_avail <- sort(unique(rel[is.finite(rel) & rel <= -2]), decreasing = TRUE)
  Lmax <- if (is.null(leads)) 5L else leads
  Kpre <- utils::head(pre_avail, Lmax)
  Kset <- c(sort(Kpre), Kpost)                    # omit -1 (reference)
  post_cols <- which(Kset %in% Kpost)

  ## Restrict to the estimation window. Treated observations at event times
  ## outside [min(Kset), max(Kset)] (leads beyond `leads`, horizons beyond the
  ## cap) must NOT silently join the comparison group -- that would contaminate
  ## the baseline. We keep never-treated units, the -1 reference period, and any
  ## event time in Kset, and drop the rest (with a note).
  in_window <- !is.finite(rel) | rel == -1L | rel %in% Kset
  n_drop <- sum(!in_window)
  if (n_drop > 0L) {
    message(sprintf(
      paste0("simpower: dropped %d obs at event times outside the window ",
             "[%d, %d] (leads > %d or horizons > %d)."),
      n_drop, min(Kset), max(Kset), if (length(Kpre)) abs(min(Kpre)) else 1L, H))
    yv  <- yv[in_window]
    uid <- as.integer(factor(uid[in_window]))
    ti  <- as.integer(factor(ti[in_window]))
    rel <- rel[in_window]                          # relative times are unchanged
    if (!is.null(W)) W <- W[in_window, , drop = FALSE]
    N <- max(uid)
  }

  ## shape weights over post horizons 0..H (peak normalised to 1)
  wpost <- switch(shape,
    constant   = rep(1, length(Kpost)),
    increasing = (Kpost + 1) / (H + 1),
    decreasing = (H + 1 - Kpost) / (H + 1)
  )
  names(wpost) <- as.character(Kpost)

  ## fixed design pieces (independent of y): event-time dummies and exposure
  rel0 <- ifelse(is.finite(rel), rel, -999L)      # never-treated -> baseline
  Dmat <- vapply(Kset, function(k) as.numeric(rel0 == k), numeric(length(rel0)))
  Dmat <- matrix(Dmat, ncol = length(Kset), dimnames = list(NULL, .knm(Kset)))
  Pexp <- ifelse(rel0 %in% Kpost, wpost[match(rel0, Kpost)], 0)
  Pexp[is.na(Pexp)] <- 0

  ## eventually-treated units = those with any finite relative event time
  treated_ids <- as.integer(sort(unique(uid[is.finite(rel)])))
  real_units  <- N
  real_treat  <- length(treated_ids)
  sd_y <- stats::sd(yv - stats::ave(yv, uid))            # within-unit SD of Y

  one_rep <- function() {
    sm <- .subsample_units(uid, treated_ids, n_units, n_treated)
    u  <- as.integer(factor(uid[sm])); t <- as.integer(factor(ti[sm]))
    ys <- yv[sm]; Ds <- Dmat[sm, , drop = FALSE]; Ps <- Pexp[sm]
    Ws <- if (is.null(W)) NULL else W[sm, , drop = FALSE]

    rec  <- .recombine_panel(ys, u, t)
    keep <- rec$keep
    yk <- rec$y[keep]
    uk <- as.integer(factor(u[keep])); tk <- as.integer(factor(t[keep]))
    Dk <- Ds[keep, , drop = FALSE]; Pk <- Ps[keep]
    Wk <- if (is.null(Ws)) NULL else Ws[keep, , drop = FALSE]

    ## overall (shaped exposure) and dynamic (per-horizon) fits share the same
    ## unit + time absorbing step, so demean once and split (see .fit_event_combined).
    fit <- .fit_event_combined(yk, Pk, Dk, uk, tk, Wk)
    bpost  <- fit$b[post_cols];  sepost <- fit$se[post_cols]

    coef <- c(overall = unname(fit$ov["b"]), stats::setNames(bpost, .knm(Kpost)))
    se   <- c(overall = unname(fit$ov["se"]), stats::setNames(sepost, .knm(Kpost)))
    ## treated units contributing a post-treatment obs after recombination
    ntr <- length(unique(uk[Pk != 0]))
    list(coef = coef, se = se,
         diag = c(units = length(unique(uk)), treated = ntr))
  }

  sim <- .run_sim(one_rep, reps, seed)

  ## overall curve
  bs <- sim$B[, "overall"]; ses <- sim$S[, "overall"]
  .warn_dropped(bs, ses, reps)
  grid  <- auto_grid(ses, emax, step)
  power <- power_from_null(bs, ses, grid, alpha, alternative, crit_method)
  mde   <- mde_targets(grid, power, power_targets)

  ## per-horizon MDE (each horizon on its own grid)
  th1 <- power_targets[1]
  hz <- lapply(seq_along(Kpost), function(j) {
    nm <- .knm(Kpost[j])
    w  <- unname(wpost[j])
    b <- sim$B[, nm]; s <- sim$S[, nm]
    g <- auto_grid(s, NULL, NULL)
    p <- power_from_null(b, s, g, alpha, alternative, crit_method)
    m <- mde_from_curve(g, p, th1)
    ## `mde` is the smallest *actual* effect detectable at this horizon on its
    ## own. Under the assumed `shape` the actual effect at horizon k is w_k * beta
    ## (beta = the peak effect), so `mde_peak = mde / w_k` is the peak effect at
    ## which this horizon becomes individually detectable -- the shape-aware
    ## per-horizon MDE, in the same peak-effect units as the overall MDE. For a
    ## "constant" shape every w_k = 1, so mde_peak == mde.
    data.frame(horizon  = Kpost[j],
               weight   = w,
               se_mean  = mean(s, na.rm = TRUE),
               mde      = m,
               mde_peak = m / w,
               mde_std  = m / sd_y)
  })
  horizon <- do.call(rbind, hz)
  attr(horizon, "target") <- th1

  req_units   <- n_units   %||% real_units
  req_treated <- n_treated %||% real_treat
  sample_check <- data.frame(
    metric  = c("units (clusters)", "treated units"),
    target  = c(req_units, req_treated),
    sim_avg = c(mean(sim$D[, "units"]), mean(sim$D[, "treated"]))
  )

  new_simpower(
    design = "event", call = cl,
    grid = grid, power = data.frame(effect = grid, power = power),
    mde = mde, bs = bs, ses = ses,
    alpha = alpha, alternative = alternative, reps = reps,
    request = list(n_units = n_units, n_treated = n_treated, controls = controls,
                   shape = shape, y = y, id = id, time = time),
    sample_check = sample_check,
    extras = list(shape = shape, horizon = horizon, H = H, Kpost = Kpost, sd_y = sd_y),
    data = data
  )
}

## Combined event-study fit. The shaped-exposure "overall" regression and the
## dynamic per-horizon regression share the same unit + time absorbing step, so
## we demean [exposure, event dummies, controls, y] once and run the two OLS
## fits on the appropriate column blocks. Numerically identical to fitting them
## separately (demeaning is linear), at roughly half the absorbing cost.
.fit_event_combined <- function(y, Pexp, Dmat, uid, ti, Xextra = NULL) {
  G <- length(unique(uid)); Tn <- length(unique(ti))
  kd <- ncol(Dmat)
  wc <- if (is.null(Xextra)) 0L else ncol(as.matrix(Xextra))
  M  <- cbind(Pexp, Dmat, Xextra, y)
  Md <- .demean_twoway(M, uid, ti)
  Pd <- Md[, 1L, drop = FALSE]
  Dd <- Md[, 1L + seq_len(kd), drop = FALSE]
  Wd <- if (wc) Md[, 1L + kd + seq_len(wc), drop = FALSE] else NULL
  yd <- Md[, ncol(Md)]

  fov <- .ols_vcov(cbind(Pd, Wd), yd, cluster = uid, dfK = G + Tn + 1L + wc)
  ov  <- if (is.null(fov)) c(b = NA_real_, se = NA_real_)
         else c(b = fov$beta[1L], se = sqrt(max(fov$V[1L, 1L], 0)))

  fdy <- .ols_vcov(cbind(Dd, Wd), yd, cluster = uid, dfK = G + Tn + kd + wc)
  if (is.null(fdy)) { b <- rep(NA_real_, kd); se <- rep(NA_real_, kd) }
  else { b <- fdy$beta[seq_len(kd)]; se <- sqrt(pmax(diag(fdy$V)[seq_len(kd)], 0)) }

  list(ov = ov, b = b, se = se)
}

## dynamic event-study fit: coefficients on the event-time dummies in Dmat,
## absorbing unit + time FE, unit-clustered SE. (Retained for direct use/tests;
## the simulation path uses .fit_event_combined.)
.fit_event_dynamic <- function(y, Dmat, uid, ti, Xextra = NULL) {
  G <- length(unique(uid)); Tn <- length(unique(ti))
  kd <- ncol(Dmat)
  M  <- cbind(Dmat, Xextra, y)
  Md <- .demean_twoway(M, uid, ti)
  p  <- kd + (if (is.null(Xextra)) 0L else ncol(as.matrix(Xextra)))
  Xd <- Md[, seq_len(p), drop = FALSE]
  yd <- Md[, p + 1L]
  fit <- .ols_vcov(Xd, yd, cluster = uid, dfK = G + Tn + p)
  if (is.null(fit)) return(list(b = rep(NA_real_, kd), se = rep(NA_real_, kd)))
  list(b = fit$beta[seq_len(kd)],
       se = sqrt(pmax(diag(fit$V)[seq_len(kd)], 0)))
}

## column names for event-time coefficients
.knm <- function(k) ifelse(k < 0, paste0("lead", abs(k)), paste0("h", k))
