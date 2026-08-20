#' Power and MDE for an instrumental-variables (2SLS) design
#'
#' Simulation-based power analysis for the coefficient on an endogenous
#' regressor `d` instrumented by `z`, estimated by two-stage least squares.
#' Because 2SLS is linear in the outcome, the analytic effect-size sweep is
#' exact: adding a true effect \eqn{\tau} shifts the 2SLS coefficient on `d` by
#' exactly \eqn{\tau} and leaves the SE unchanged. Instrument strength (the
#' first stage) is preserved because the real `z` and `d` are used in every
#' simulated fit --- weak instruments show up as low power automatically.
#'
#' @inheritParams power_twfe
#' @param d Column name of the endogenous regressor whose coefficient is tested.
#' @param z Column name(s) of the excluded instrument(s).
#' @param id Optional column name of a clustering / panel unit id. If supplied
#'   (with `time`) the outcome is recombined by donor series to preserve
#'   within-unit autocorrelation and SEs are clustered by `id`; otherwise a
#'   cross-sectional residual-permutation null with HC1 SEs is used.
#' @param time Optional time column (only used, with `id`, for panel-style
#'   recombination).
#' @param n_units,n_treated Number of units and treated units (`d == 1`) to draw
#'   per simulated sample. Units are clusters when `id` is given, otherwise
#'   individual observations.
#'
#' @param L Optional column name switching on **hypothetical-instrument mode**
#'   (supply either `z` or `L`, not both). `L` is any variable whose
#'   distributional character you want the imagined instrument to inherit. On
#'   every simulation rep the `L` series is randomly permuted -- breaking its
#'   correlation with both the treatment and the outcome -- and a synthetic
#'   instrument is built as \eqn{z_r = r\,\tilde X + \sqrt{1-r^2}\,\tilde
#'   L_\perp}, where \eqn{\tilde X} is the standardised treatment and
#'   \eqn{\tilde L_\perp} the permuted `L` residualised on the treatment and
#'   standardised. By construction the *in-sample Pearson correlation between
#'   the instrument and the treatment is exactly* `r` on every rep, so the
#'   resulting curves isolate what instrument strength alone does to power.
#'   One power curve and one MDE row are returned per value of `rho`.
#' @param rho Numeric vector of target first-stage correlations (Pearson r
#'   between the synthetic instrument and the treatment), each strictly
#'   between 0 and 1. Only used in hypothetical-instrument mode. Default
#'   `c(0.1, 0.3, 0.5)`.
#'
#' @return An object of class `"simpower"` (design `"iv"`). `extras` holds the
#'   real-data 2SLS estimate and the first-stage F statistic. In
#'   hypothetical-instrument mode `$power` holds one curve per `rho` (long
#'   format, with a `rho` column), `$mde` is a data frame with one row per
#'   `rho`, and `extras$first_stage_F` reports the implied first-stage F at
#'   each `rho` (deterministic given the sample size, because the in-sample
#'   correlation is exact by construction).
#'
#' @examples
#' set.seed(1)
#' df <- sim_iv(n = 800, strength = 0.4)
#' pw <- power_iv(df, y = "y", d = "d", z = "z", reps = 200)
#' pw
#'
#' ## hypothetical instrument: how would power look with cor(z, d) = .1/.3/.5?
#' df$candidate <- rnorm(nrow(df))
#' ph <- power_iv(df, y = "y", d = "d", L = "candidate",
#'                rho = c(0.1, 0.3, 0.5), reps = 200)
#' ph
#' @export
power_iv <- function(data, y, d, z = NULL, id = NULL, time = NULL, controls = NULL,
                     n_units = NULL, n_treated = NULL,
                     L = NULL, rho = c(0.1, 0.3, 0.5),
                     reps = 1000, emax = NULL, step = NULL,
                     alpha = 0.05, alternative = c("greater", "two.sided", "less"),
                     crit_method = c("empirical", "normal"),
                     power_targets = c(0.8, 0.9), seed = 5000, het = NULL) {
  alternative <- match.arg(alternative)
  crit_method <- match.arg(crit_method)
  cl <- match.call()
  hspec <- .het_check(het)
  henv  <- if (!is.null(hspec)) .het_new_stream(seed)

  hypo <- !is.null(L)
  if (hypo && !is.null(z)) {
    stop("Supply either `z` (a real instrument) or `L` (hypothetical-instrument mode), not both.",
         call. = FALSE)
  }
  if (!hypo && is.null(z)) {
    stop("Supply `z` (instrument column name(s)), or `L` for hypothetical-instrument mode.",
         call. = FALSE)
  }
  if (hypo) {
    rho <- sort(unique(as.numeric(rho)))
    if (!length(rho) || any(!is.finite(rho)) || any(rho <= 0) || any(rho >= 1)) {
      stop("`rho` must contain correlations strictly between 0 and 1.", call. = FALSE)
    }
    if (!is.null(time) && !is.null(id)) {
      stop(paste0("Hypothetical-instrument mode is cross-sectional: omit `time` ",
                  "(`id` may still be supplied for clustered SEs)."), call. = FALSE)
    }
  }

  yv <- .numcol(data, y, "y")
  dv <- .numcol(data, d, "d")
  Z  <- if (hypo) NULL else {
    Zm <- vapply(z, function(zz) .numcol(data, zz, "z"), numeric(nrow(data)))
    matrix(Zm, ncol = length(z), dimnames = list(NULL, z))
  }
  Lv <- if (hypo) .numcol(data, L, "L") else NULL
  W  <- .controls_matrix(data, controls)
  has_panel <- !is.null(id) && !is.null(time)
  idv <- if (!is.null(id)) .col(data, id, "id") else NULL
  tvv <- if (!is.null(time)) .col(data, time, "time") else NULL

  keep0 <- stats::complete.cases(yv, dv,
                                 if (is.null(Z)) NULL else Z,
                                 if (is.null(Lv)) NULL else Lv,
                                 if (is.null(W)) NULL else W,
                                 if (is.null(idv)) NULL else idv,
                                 if (is.null(tvv)) NULL else tvv)
  yv <- yv[keep0]; dv <- dv[keep0]
  if (!is.null(Z))  Z  <- Z[keep0, , drop = FALSE]
  if (!is.null(Lv)) Lv <- Lv[keep0]
  if (!is.null(W)) W <- W[keep0, , drop = FALSE]

  n <- length(yv)
  uid <- if (!is.null(idv)) as.integer(factor(idv[keep0])) else seq_len(n)
  ti  <- if (!is.null(tvv)) as.integer(factor(tvv[keep0])) else NULL

  ## treated units. For a binary endogenous regressor these are units with
  ## d == 1; for a continuous one ("dose") there is no natural 0/1 split, so we
  ## label units with above-median average exposure as "treated" purely so that
  ## `n_treated` subsampling and the diagnostic remain meaningful.
  is_bin <- all(dv %in% c(0, 1))
  med <- stats::median(dv)
  is_trt_unit <- function(v) if (is_bin) any(v == 1) else mean(v) > med
  trt_flag <- tapply(dv, uid, is_trt_unit)
  treated_ids <- as.integer(names(trt_flag))[trt_flag %in% TRUE]
  real_units <- length(unique(uid)); real_treat <- length(treated_ids)
  sd_y <- if (has_panel) stats::sd(yv - stats::ave(yv, uid)) else stats::sd(yv)

  ## exogenous block for the cross-sectional residual null
  Eexog <- W

  one_rep <- function() {
    ss <- .subsample_units(uid, treated_ids, n_units, n_treated)
    sm <- ss$rows
    ys <- yv[sm]; ds <- dv[sm]
    Zs <- if (is.null(Z)) NULL else Z[sm, , drop = FALSE]
    Ws <- if (is.null(W)) NULL else W[sm, , drop = FALSE]
    us <- as.integer(factor(ss$uid))

    if (hypo) {
      ## Hypothetical-instrument mode. Impose the null on y (residual
      ## permutation, as in the cross-sectional path), permute L to break its
      ## link to both x and y, then build one exact-correlation synthetic
      ## instrument per requested rho and fit 2SLS with each. The same
      ## recombined y and permuted L are shared across the rho values, so the
      ## curves differ only through instrument strength (common random numbers).
      Ek  <- if (is.null(Ws)) NULL else Ws
      yk  <- .recombine_resid(ys, Ek)
      Lp  <- sample(Lv[sm])
      clu <- if (!is.null(idv)) us else NULL
      inj <- if (is.null(hspec)) NULL else
        .het_draw(hspec, length(unique(us)), henv)[us] * ds
      tmpl <- if (is.null(hspec)) c(b = 0, se = 0)
              else c(b = 0, se = 0, lam = 0, v1 = 0, v2 = 0)
      est <- vapply(rho, function(r) {
        zr <- .synth_instrument(ds, Lp, r)
        if (is.null(zr)) return(tmpl + NA_real_)
        .fit_tsls(yk, ds, matrix(zr, ncol = 1L), Ws, cluster = clu, inject = inj)
      }, tmpl)
      nm <- .rnm(rho)
      return(list(coef = stats::setNames(est["b", ], nm),
                  se   = stats::setNames(est["se", ], nm),
                  diag = c(units = length(unique(us)),
                           treated = sum(tapply(ds, us, is_trt_unit)),
                           nobs = length(ys)),
                  het  = if (!is.null(hspec)) stats::setNames(
                    as.vector(est[c("lam", "v1", "v2"), ]),
                    paste(rep(nm, each = 3L), c("lam", "v1", "v2"), sep = "."))))
    }

    if (has_panel) {
      ts  <- as.integer(factor(ti[sm]))
      rec <- .recombine_panel_ctrl(ys, us, ts, Ws, fe = FALSE)  # .fit_tsls has no FEs
      keep <- rec$keep
      yk <- rec$y[keep]; dk <- ds[keep]; Zk <- Zs[keep, , drop = FALSE]
      Wk <- if (is.null(Ws)) NULL else Ws[keep, , drop = FALSE]
      uk <- as.integer(factor(us[keep]))
      inj <- if (is.null(hspec)) NULL else
        .het_draw(hspec, length(unique(uk)), henv)[uk] * dk
      est <- .fit_tsls(yk, dk, Zk, Wk, cluster = uk, inject = inj)
      ntreat <- sum(tapply(dk, uk, is_trt_unit))
      nunit  <- length(unique(uk))
    } else {
      Ek <- if (is.null(Ws)) NULL else Ws
      yk <- .recombine_resid(ys, Ek)
      clu <- if (!is.null(idv)) us else NULL
      inj <- if (is.null(hspec)) NULL else
        .het_draw(hspec, length(unique(us)), henv)[us] * ds
      est <- .fit_tsls(yk, ds, Zs, Ws, cluster = clu, inject = inj)
      ntreat <- sum(tapply(ds, us, is_trt_unit))
      nunit  <- length(unique(us))
    }
    list(coef = c(effect = unname(est["b"])),
         se   = c(effect = unname(est["se"])),
         diag = c(units = nunit, treated = ntreat),
         het  = if (!is.null(hspec)) est[c("lam", "v1", "v2")])
  }

  sim <- .run_sim(one_rep, reps, seed)

  sample_check <- data.frame(
    metric  = c(if (!is.null(id)) "units (clusters)" else "observations", "treated units"),
    target  = c(n_units %||% real_units, n_treated %||% real_treat),
    sim_avg = c(mean(sim$D[, "units"]), mean(sim$D[, "treated"]))
  )

  if (hypo) {
    ## One curve + MDE row per rho, on a common effect grid (scaled to the
    ## weakest instrument so every curve is fully visible).
    .warn_dropped(sim$B, sim$S, reps * length(rho))
    grid <- auto_grid(sim$S[, 1L], emax, step, widen = !is.null(hspec))
    curves <- lapply(seq_along(rho), function(j) {
      pj <- if (is.null(hspec)) {
        power_from_null(sim$B[, j], sim$S[, j], grid,
                        alpha, alternative, crit_method)
      } else {
        power_from_null_het(sim$B[, j], sim$S[, j],
                            sim$H[, 3L * (j - 1L) + 1:3, drop = FALSE], grid,
                            alpha, alternative, crit_method)
      }
      data.frame(effect = grid, power = pj, rho = rho[j])
    })
    power_df <- do.call(rbind, curves)
    mde_df <- cbind(
      data.frame(rho = rho),
      do.call(rbind, lapply(seq_along(rho), function(j) {
        as.data.frame(as.list(mde_targets(grid, curves[[j]]$power, power_targets)),
                      check.names = FALSE)
      }))
    )
    ## Implied first-stage F at each rho: the in-sample cor(z, x) is exact, so
    ## for a single instrument F = t^2 = r^2 (n - 2) / (1 - r^2).
    nbar <- mean(sim$D[, "nobs"])
    Fimp <- stats::setNames(rho^2 * (nbar - 2) / (1 - rho^2), .rnm(rho))

    return(new_simpower(
      design = "iv", call = cl,
      grid = grid, power = power_df,
      mde = mde_df, bs = sim$B, ses = sim$S,
      alpha = alpha, alternative = alternative, reps = reps,
      request = list(n_units = n_units, n_treated = n_treated, controls = controls,
                     y = y, d = d, L = L, rho = rho, id = id),
      sample_check = sample_check,
      extras = list(estimate = c(b = NA_real_, se = NA_real_),
                    first_stage_F = Fimp, rho = rho, L = L,
                    clustered = !is.null(idv), panel = FALSE, sd_y = sd_y,
                    het = if (!is.null(hspec)) list(label = .het_label(hspec),
                                                    spec = hspec$spec)),
      data = data
    ))
  }

  bs <- sim$B[, 1]; ses <- sim$S[, 1]
  .warn_dropped(bs, ses, reps)
  grid  <- auto_grid(ses, emax, step, widen = !is.null(hspec))
  hc    <- .het_curves(bs, ses, sim$H, grid, alpha, alternative, crit_method,
                       power_targets, hspec, auto_extend = is.null(emax))
  power <- hc$power
  mde   <- hc$mde
  grid  <- hc$grid

  ## real-data reference: 2SLS estimate + first-stage F (diagnostics only)
  ref <- .iv_reference(yv, dv, Z, W, if (!is.null(idv)) uid else NULL)

  new_simpower(
    design = "iv", call = cl,
    grid = grid, power = data.frame(effect = grid, power = power),
    mde = mde, bs = bs, ses = ses,
    alpha = alpha, alternative = alternative, reps = reps,
    request = list(n_units = n_units, n_treated = n_treated, controls = controls,
                   y = y, d = d, z = z, id = id, time = time),
    sample_check = sample_check,
    extras = list(estimate = ref$estimate, first_stage_F = ref$F,
                  clustered = !is.null(idv), panel = has_panel, sd_y = sd_y,
                  het = hc$extras),
    data = data
  )
}

## Synthetic instrument with *exact* in-sample correlation r to the treatment x:
## z = r * x_std + sqrt(1 - r^2) * e_std, where e is the permuted carrier
## series `Lp` residualised on (1, x) -- exactly orthogonal to x in-sample --
## and both pieces are standardised. Then cor(z, x) = r identically. Returns
## NULL when x or the residualised carrier is degenerate (constant), so the
## caller can drop the rep.
.synth_instrument <- function(x, Lp, r) {
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= 0) return(NULL)
  xs <- (x - mean(x)) / sx
  e  <- stats::lm.fit(cbind(1, xs), Lp)$residuals
  se <- stats::sd(e)
  if (!is.finite(se) || se <= 0) return(NULL)
  r * xs + sqrt(1 - r^2) * (e / se)
}

## column names for the per-rho coefficients ("r0.1", "r0.3", ...)
.rnm <- function(rho) paste0("r", format(rho, trim = TRUE, drop0trailing = TRUE))

## Real-data 2SLS point estimate (via AER if available) and first-stage F.
.iv_reference <- function(y, d, Z, W, cluster = NULL) {
  have_pkgs <- requireNamespace("AER", quietly = TRUE) &&
               requireNamespace("sandwich", quietly = TRUE)
  est <- if (!have_pkgs) c(b = NA_real_, se = NA_real_) else tryCatch({
    df <- data.frame(y = y, d = d)
    for (j in seq_len(ncol(Z))) df[[paste0("z", j)]] <- Z[, j]
    znames <- paste0("z", seq_len(ncol(Z)))
    wnames <- character(0)
    if (!is.null(W)) { for (j in seq_len(ncol(W))) { nm <- paste0("w", j); df[[nm]] <- W[, j]; wnames <- c(wnames, nm) } }
    rhs  <- paste(c("d", wnames), collapse = " + ")
    inst <- paste(c(znames, wnames), collapse = " + ")
    fml  <- stats::as.formula(sprintf("y ~ %s | %s", rhs, inst))
    m <- AER::ivreg(fml, data = df)
    vc <- if (!is.null(cluster)) sandwich::vcovCL(m, cluster = cluster) else sandwich::vcovHC(m, type = "HC1")
    c(b = unname(stats::coef(m)["d"]), se = sqrt(vc["d", "d"]))
  }, error = function(e) c(b = NA_real_, se = NA_real_))

  Fstat <- tryCatch({
    ones <- rep(1, length(d))
    Xfull <- cbind(ones, W, Z); Xrest <- cbind(ones, W)
    fs_full <- stats::lm.fit(Xfull, d)
    fs_rest <- stats::lm.fit(Xrest, d)
    rss1 <- sum(fs_full$residuals^2); rss0 <- sum(fs_rest$residuals^2)
    q <- ncol(Z); nobs <- length(d); kk <- ncol(Xfull)
    ((rss0 - rss1) / q) / (rss1 / (nobs - kk))
  }, error = function(e) NA_real_)

  list(estimate = est, F = Fstat)
}
