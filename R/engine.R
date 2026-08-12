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
## Draw a design with `n_units` total and `n_treated` treated units. Returns a
## logical row mask over the working data. When both are NULL, keep everything.
.subsample_units <- function(uid, treated_ids, n_units = NULL, n_treated = NULL) {
  ## Common default: no subsampling requested -> keep everything. Return before
  ## doing any pool bookkeeping, since this runs once per simulation rep.
  if (is.null(n_units) && is.null(n_treated)) {
    return(rep(TRUE, length(uid)))
  }
  units   <- unique(uid)
  is_trt  <- units %in% treated_ids
  trt_pool <- units[is_trt]
  ctl_pool <- units[!is_trt]

  n_units <- n_units %||% length(units)
  if (is.null(n_treated)) {
    # keep the natural treated share among the sampled units
    share <- length(trt_pool) / length(units)
    n_treated <- round(share * n_units)
  }
  n_treated <- min(n_treated, length(trt_pool), n_units)
  n_ctl     <- min(n_units - n_treated, length(ctl_pool))
  keep_units <- c(
    if (n_treated > 0) sample(trt_pool, n_treated) else integer(0),
    if (n_ctl     > 0) sample(ctl_pool, n_ctl)     else integer(0)
  )
  uid %in% keep_units
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
  B  <- matrix(NA_real_, reps, m)
  S  <- matrix(NA_real_, reps, m)
  D  <- matrix(NA_real_, reps, max(nd, 1L))
  B[1L, ] <- first$coef; S[1L, ] <- first$se
  if (nd) D[1L, ] <- first$diag
  if (reps >= 2L) {
    for (r in 2:reps) {
      est <- one_rep()
      B[r, ] <- est$coef; S[r, ] <- est$se
      if (nd) D[r, ] <- est$diag
    }
  }
  colnames(B) <- names(first$coef)
  colnames(S) <- names(first$coef)
  if (nd) colnames(D) <- names(first$diag)
  list(B = B, S = S, D = D)
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
