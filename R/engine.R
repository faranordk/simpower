## ---------------------------------------------------------------------------
## Recombination (imposing the empirical null) and unit subsampling.
## ---------------------------------------------------------------------------

## Panel recombination: give each unit a *different* unit's outcome series,
## keeping the time alignment and each unit's own regressors. This destroys any
## true relationship while preserving each series' autocorrelation, variance,
## and the panel's clustering/missingness structure. Rows whose donor cell is
## missing (unbalanced panels) are dropped, so a realistic missingness pattern
## is retained. Returns a logical `keep` mask and the recombined outcome.
.recombine_panel <- function(y, uid, ti) {
  N <- max(uid); Tn <- max(ti)
  Ymat <- matrix(NA_real_, N, Tn)
  Ymat[cbind(uid, ti)] <- y
  donor <- .derangement(N)[uid]          # a random *other* donor unit per unit
  yp <- Ymat[cbind(donor, ti)]
  keep <- !is.na(yp)
  list(y = yp, keep = keep)
}

## Control-preserving panel recombination.
##
## `.recombine_panel()` swaps each unit's *entire* outcome series for a donor's,
## which destroys not only the y ~ treatment link (intended) but also the
## y ~ controls link (unintended). Estimators that partial out controls then
## face a recombined outcome the controls cannot predict, so the simulated
## residual variance -- and with it the SE and the MDE -- stays at its
## no-controls level: the efficiency gain from controls that predict Y is lost.
##
## Fix: estimate the control-predicted component of y under the two-way FE
## structure (a within-regression of y on the controls only -- the tested
## regressor is deliberately excluded, so any true-effect variance stays in the
## residual and is scrambled into noise, per the empirical-null philosophy),
## recombine the *control-purged* residual series across units (preserving its
## autocorrelation, variance, and missingness), and add each row's OWN
## control-predicted component back. The refitted controls now explain exactly
## the share of variance they explain in the real data, while the tested
## regressor still faces a pure-noise outcome. With `W = NULL` this is
## byte-identical to `.recombine_panel()` (same RNG stream: the partialling
## draws no random numbers).
##
## `fe` must match the refitting estimator's own specification: TRUE (the
## default) estimates gamma-hat under absorbed unit + time fixed effects, as
## the TWFE / event-study estimators do; FALSE estimates it by a pooled
## regression of y on (1, W), matching estimators that include no fixed
## effects (the panel-IV `.fit_tsls()`). Using the wrong block would leave a
## W-predictable component in the recombined residual (or remove one the
## estimator cannot see), distorting the simulated residual variance.
.recombine_panel_ctrl <- function(y, uid, ti, W = NULL, fe = TRUE) {
  cf <- .panel_control_fit(y, W, uid, ti, fe = fe)
  if (is.null(cf)) return(.recombine_panel(y, uid, ti))
  rec <- .recombine_panel(y - cf, uid, ti)
  rec$y <- rec$y + cf                     # rows with a missing donor stay NA
  rec
}

## The control-predicted component of y, on the raw scale (W %*% gamma-hat,
## no tested regressor). With `fe = TRUE`, gamma-hat comes from the two-way
## within regression of y on W; with `fe = FALSE`, from the pooled regression
## of y on (1, W) (the intercept is dropped from the returned component -- a
## constant shift is irrelevant to recombination and estimation alike). NULL
## when there are no controls or the cross-product is singular (callers then
## fall back to the uncorrected recombination).
.panel_control_fit <- function(y, W, uid, ti, fe = TRUE) {
  if (is.null(W)) return(NULL)
  W <- as.matrix(W)
  k <- ncol(W)
  if (fe) {
    Md <- .demean_twoway(cbind(W, y), uid, ti)
    Wd <- Md[, seq_len(k), drop = FALSE]
    yd <- Md[, k + 1L]
  } else {
    Wd <- sweep(W, 2L, colMeans(W))       # pooled: demean = intercept absorbed
    yd <- y - mean(y)
  }
  Ainv <- .safe_inv(crossprod(Wd))
  if (is.null(Ainv)) return(NULL)
  as.vector(W %*% (Ainv %*% crossprod(Wd, yd)))
}

## A random derangement of 1..N (no element maps to itself), so no unit is ever
## paired with its own outcome series. A plain permutation leaves ~1 fixed point
## on average, which would leak a little true signal into the imposed null.
## Uses Sattolo's algorithm (an O(N) draw of a single-cycle permutation, which is
## always a derangement for N >= 2).
.derangement <- function(N) {
  if (N < 2L) return(seq_len(N))
  p <- seq_len(N)
  for (i in N:2L) {
    j <- sample.int(i - 1L, 1L)          # j strictly less than i
    tmp <- p[i]; p[i] <- p[j]; p[j] <- tmp
  }
  p
}

## Cross-sectional recombination: preserve the part of y explained by the
## exogenous columns `Eexog` (controls, running-variable trend, intercept) and
## permute the residual. This imposes the null on the *tested* regressor while
## keeping the exogenous structure and the residual variance intact. If Eexog
## is NULL, this reduces to a plain permutation of y.
.recombine_resid <- function(y, Eexog = NULL) {
  if (is.null(Eexog)) return(sample(y))
  E <- cbind(1, as.matrix(Eexog))
  Ainv <- .safe_inv(crossprod(E))
  if (is.null(Ainv)) return(sample(y))
  fit  <- E %*% (Ainv %*% crossprod(E, y))
  res  <- y - fit
  as.vector(fit + sample(res))
}

## Classify units as treated / control at the unit level for subsampling.
## `treated_ids` is a vector of unit ids considered treated.
## Draw a design with `n_units` total and `n_treated` treated units. Returns
## `list(rows, uid)`: row indices into the working data (with repeated blocks
## when upsampling) and a matching unit-id vector in which each drawn unit
## *instance* is a distinct unit. When both are NULL, keep everything.
##
## Either count may EXCEED what the data contain ("would collecting more data
## help?"). The surplus is drawn by replicating units -- each pool enters
## floor(n/pool) times in full plus a random without-replacement draw of the
## remainder -- and every copy becomes its own cluster (a cluster bootstrap).
## Copies keep their complete outcome series and regressors, so the
## extrapolation assumes additional units would look like those already
## observed. (In the panel recombination a copy can occasionally be assigned
## its twin's -- i.e. its own -- outcome series; with N units this touches
## O(1/N) of the sample and is negligible.) When the request is within the
## sample, the drawn rows and the RNG stream are identical to previous
## versions, so existing results reproduce exactly.
.subsample_units <- function(uid, treated_ids, n_units = NULL, n_treated = NULL) {
  ## Common default: no subsampling requested -> keep everything. Return before
  ## doing any pool bookkeeping, since this runs once per simulation rep.
  if (is.null(n_units) && is.null(n_treated)) {
    return(list(rows = seq_along(uid), uid = uid))
  }
  units   <- unique(uid)
  is_trt  <- units %in% treated_ids
  trt_pool <- units[is_trt]
  ctl_pool <- units[!is_trt]

  if (is.null(n_units)) {
    ## `n_treated` alone beyond the treated pool: keep every control unit and
    ## grow the total, rather than silently crowding the controls out.
    n_units <- if (is.null(n_treated)) length(units)
               else max(length(units), n_treated + length(ctl_pool))
  }
  if (is.null(n_treated)) {
    # keep the natural treated share among the sampled units
    share <- length(trt_pool) / length(units)
    n_treated <- round(share * n_units)
  }
  n_treated <- min(n_treated, n_units)
  if (!length(trt_pool)) n_treated <- 0L
  n_ctl <- n_units - n_treated
  if (!length(ctl_pool)) n_ctl <- 0L
  draws <- c(.draw_units(trt_pool, n_treated), .draw_units(ctl_pool, n_ctl))

  if (!anyDuplicated(draws)) {
    ## within-sample draw: same rows (ascending) as the old mask-based code
    rows <- which(uid %in% draws)
    return(list(rows = rows, uid = uid[rows]))
  }
  ## upsampled draw: each instance becomes its own unit (cluster bootstrap)
  by_unit <- split(seq_along(uid), uid)
  rl <- by_unit[as.character(draws)]
  list(rows = unlist(rl, use.names = FALSE),
       uid  = rep(seq_along(rl), lengths(rl)))
}

## Draw `n` unit ids from `pool`: a plain subsample when n <= |pool| (same RNG
## usage as sample()), otherwise floor(n/|pool|) full copies of the pool plus a
## without-replacement draw of the remainder.
.draw_units <- function(pool, n) {
  m <- length(pool)
  if (n <= 0L || !m) return(pool[0])
  if (n <= m) return(pool[sample.int(m, n)])
  k <- n %/% m; r <- n %% m
  c(rep(pool, k), if (r > 0L) pool[sample.int(m, r)] else pool[0])
}

## ---------------------------------------------------------------------------
## Heterogeneous treatment effects (the `het` option).
##
## The default sweep injects the SAME effect e into every unit. With `het`,
## each unit instead receives e * m_i, where the multipliers m_i are drawn
## per simulation rep from a user-chosen mean-one distribution -- so e remains
## the AVERAGE effect and the power curves are directly comparable to the
## constant-effect ones. Because every estimator here is linear in y and every
## SE is a sandwich, no refitting across effect sizes is needed: injecting
## e * d (d = m * treatment direction) shifts the tested coefficient by
## exactly e * lambda, with lambda the estimator applied to d, and turns the
## sandwich variance into an exact quadratic v0 + e*v1 + e^2*v2 (residuals are
## linear in e, the meat is quadratic). Each fit therefore returns three extra
## scalars per rep and the whole heterogeneous power curve stays a cheap
## post-processing sweep -- exact, not bootstrapped.
## ---------------------------------------------------------------------------

## Validate a `het` spec once, up front. NULL -> NULL (constant effects).
## A numeric vector is an equal-probability multiset of effect multipliers,
## normalised to mean one (e.g. c(0, 2): the effect is absent in half the
## units and doubled in the other half). A function(n) must return n
## multipliers with mean ~ 1 (checked on a large draw).
.het_check <- function(het) {
  if (is.null(het)) return(NULL)
  if (is.numeric(het)) {
    if (!length(het) || any(!is.finite(het))) {
      stop("`het` must be a finite numeric vector of effect multipliers.", call. = FALSE)
    }
    if (mean(het) <= 0) {
      stop("`het` multipliers must have a positive mean (e is the average effect).",
           call. = FALSE)
    }
    v <- het / mean(het)
    return(list(type = "multiset", values = v, spec = het,
                constant = length(unique(v)) == 1L))
  }
  if (is.function(het)) {
    probe <- het(1000L)
    if (!is.numeric(probe) || length(probe) != 1000L) {
      stop("`het` function must return n numeric multipliers.", call. = FALSE)
    }
    if (abs(mean(probe) - 1) > 0.05) {
      warning(sprintf(paste0("`het` function draws have mean %.3f; e is interpreted ",
                             "as the average effect, so the multipliers should have ",
                             "mean 1."), mean(probe)), call. = FALSE)
    }
    return(list(type = "function", fn = het, spec = het, constant = FALSE))
  }
  stop("`het` must be NULL, a numeric vector of multipliers, or a function(n).",
       call. = FALSE)
}

## Draw n per-unit multipliers from a checked spec, on an ISOLATED RNG stream
## (`env` holds it) so the recombination draws -- and therefore the simulated
## null estimates and SEs -- are bit-identical to a run without `het`. The
## power curves then differ only through the injected heterogeneity (common
## random numbers). A degenerate spec (all multipliers equal) draws nothing.
.het_new_stream <- function(seed) {
  env <- new.env(parent = emptyenv())
  main <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  set.seed(if (is.null(seed)) 104729L else seed + 104729L)
  env$state <- get(".Random.seed", globalenv())
  if (!is.null(main)) assign(".Random.seed", main, globalenv())
  env
}
.het_draw <- function(spec, n, env) {
  if (spec$type == "multiset" && spec$constant) return(rep(1, n))
  main <- get(".Random.seed", globalenv())
  assign(".Random.seed", env$state, globalenv())
  m <- if (spec$type == "multiset") {
    spec$values[sample.int(length(spec$values), n, replace = TRUE)]
  } else {
    spec$fn(n)
  }
  env$state <- get(".Random.seed", globalenv())
  assign(".Random.seed", main, globalenv())
  m
}

## Shared simulation runner.
##
## `one_rep()` is a design-specific closure that draws one recombined sample
## and returns a list with:
##   coef : named numeric vector of tested coefficient(s) under the null,
##   se   : matching standard error(s),
##   diag : named numeric vector of sample-composition diagnostics.
## For scalar designs `coef`/`se` are length-1; for the event study they are
## one entry per horizon.
.run_sim <- function(one_rep, reps, seed) {
  if (!is.null(seed)) set.seed(seed)
  first <- one_rep()
  m  <- length(first$coef)
  nd <- length(first$diag)
  nh <- length(first$het)                 # optional heterogeneity pieces
  B  <- matrix(NA_real_, reps, m)
  S  <- matrix(NA_real_, reps, m)
  D  <- matrix(NA_real_, reps, max(nd, 1L))
  H  <- if (nh) matrix(NA_real_, reps, nh) else NULL
  B[1L, ] <- first$coef; S[1L, ] <- first$se
  if (nd) D[1L, ] <- first$diag
  if (nh) H[1L, ] <- first$het
  if (reps >= 2L) {
    for (r in 2:reps) {
      est <- one_rep()
      B[r, ] <- est$coef; S[r, ] <- est$se
      if (nd) D[r, ] <- est$diag
      if (nh) H[r, ] <- est$het
    }
  }
  colnames(B) <- names(first$coef)
  colnames(S) <- names(first$coef)
  if (nd) colnames(D) <- names(first$diag)
  if (nh) colnames(H) <- names(first$het)
  list(B = B, S = S, D = D, H = H)
}

## Warn once if a non-trivial share of simulation reps produced an unusable
## (non-finite or zero-SE) estimate and were dropped from the power calculation
## -- e.g. because a subsample was too small to identify the design. Silent
## dropping would make the curve look more reliable than it is.
.warn_dropped <- function(bs, ses, reps, threshold = 0.05) {
  bad <- !(is.finite(bs) & is.finite(ses) & ses > 0)
  frac <- mean(bad)
  if (frac > threshold) {
    warning(sprintf(
      paste0("%.0f%% of simulation reps (%d/%d) produced an unusable estimate ",
             "and were dropped (often too few units/clusters or degrees of ",
             "freedom). The power curve uses the remaining reps."),
      100 * frac, sum(bad), length(bs)), call. = FALSE)
  }
  invisible(frac)
}

## Construct the returned S3 object.
new_simpower <- function(design, call, grid, power, mde, bs, ses,
                         alpha, alternative, reps, request,
                         sample_check, extras = list(), data = NULL) {
  structure(
    list(
      design       = design,
      call         = call,
      data         = data,      # the fitting data, so plan_mde() can re-run safely
      grid         = grid,
      power        = power,
      mde          = mde,
      null         = list(b_mean = mean(bs, na.rm = TRUE),
                          se_mean = mean(ses, na.rm = TRUE),
                          bs = bs, ses = ses),
      alpha        = alpha,
      alternative  = alternative,
      reps         = reps,
      request      = request,
      sample_check = sample_check,
      extras       = extras
    ),
    class = "simpower"
  )
}

## Small helper: pull a column by name from a data.frame with a clear error.
.col <- function(data, name, what) {
  if (is.null(name)) stop(sprintf("`%s` must be supplied.", what), call. = FALSE)
  if (!name %in% names(data)) {
    stop(sprintf("Column '%s' (%s) not found in `data`.", name, what), call. = FALSE)
  }
  data[[name]]
}

## Pull a column that must be numeric-valued (outcome, treatment, instrument,
## running variable, ...). Unlike `as.numeric(.col(...))`, this refuses to
## silently turn a factor/character into its integer level codes -- a common
## source of quietly-wrong results.
.numcol <- function(data, name, what) {
  v <- .col(data, name, what)
  if (is.factor(v) || is.character(v)) {
    stop(sprintf(paste0("Column '%s' (%s) must be numeric, but is %s. ",
                        "Convert it to a numeric coding before calling."),
                 name, what, class(v)[1L]), call. = FALSE)
  }
  as.numeric(v)
}

## Build a numeric control matrix from a character vector of column names.
## Factors/characters are expanded via model.matrix (dummy coding).
.controls_matrix <- function(data, controls) {
  if (is.null(controls) || !length(controls)) return(NULL)
  miss <- setdiff(controls, names(data))
  if (length(miss)) {
    stop(sprintf("Control column(s) not found: %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  fml <- stats::reformulate(controls)
  ## Keep every row: the default na.action drops rows with missing controls,
  ## which would misalign the returned matrix against `data`. Build the model
  ## frame with na.pass so NA rows survive as NA (removed later by the caller's
  ## complete.cases() step). model.matrix() alone does not honour na.action here.
  mf  <- stats::model.frame(fml, data = data, na.action = stats::na.pass)
  mm  <- stats::model.matrix(fml, data = mf)
  mm  <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (!ncol(mm)) return(NULL)
  mm
}
