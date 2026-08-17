# simpower 0.2.1

* **`n_units` / `n_treated` may now exceed the sample ("would collecting more
  data help?").** Previously a request beyond the observed number of units was
  silently capped at what the data contain, so sample-size planning could only
  sweep downward. Requests beyond the pool are now met by replicating observed
  units -- each pool enters `floor(n/pool)` times in full plus a random draw
  of the remainder -- with every copy treated as its own cluster (a cluster
  bootstrap). This works in all four designs and, since `plan_mde()` simply
  re-runs the fitted call, `plan_mde(pw, n_units = c(48, 96, 192))` now traces
  the MDE into hypothetical larger samples. The extrapolation assumes new
  units resemble those already observed. Within-sample requests draw exactly
  the same units (and consume the same RNG stream) as before, so existing
  results reproduce bit-for-bit. Supplying only `n_treated` beyond the treated
  pool now grows the total so all control units are kept, instead of crowding
  controls out.

* **Panel designs now capture the efficiency gain from controls that predict
  the outcome.** The panel recombination (`power_twfe()`, `power_event()`, and
  `power_iv()` with `id` + `time`) swapped each unit's entire outcome series
  for a donor's, which broke not only the y ~ treatment link (intended) but
  also the y ~ controls link (unintended): in every simulated fit the controls
  explained ~nothing, so the simulated SEs -- and the MDE -- ignored the
  variance reduction that predictive controls deliver in the real fit. The
  recombination now partials the controls' predicted component out of the
  outcome, recombines the control-purged residual series across units, and
  adds each row's own control-predicted component back. The component is
  estimated to match each estimator's own specification: a two-way within
  regression of y on the controls for `power_twfe()` / `power_event()` (whose
  estimators absorb both fixed effects), and a pooled regression of y on the
  controls for the panel branch of `power_iv()` (whose 2SLS estimator includes
  no fixed effects). Controls therefore explain exactly the
  share of variance in the simulated fits that they explain in the real data;
  the null estimates stay centered at zero and their spread shrinks in step
  with the SEs, so the null t-statistics -- and the test's size -- are
  unchanged. With `controls = NULL` the code path and RNG stream are
  identical, so results without controls are exactly reproduced. The
  cross-sectional designs (`power_rdd()`, cross-sectional `power_iv()`)
  already preserved the controls' fit via the residual-permutation null and
  are unchanged.

# simpower 0.2.0

* **New: hypothetical-instrument mode in `power_iv()`.** Supply `L =` (a column
  whose distributional character the imagined instrument inherits) instead of
  `z =`, plus a vector of target first-stage correlations
  `rho = c(0.1, 0.3, 0.5)`. Each rep permutes the `L` series (breaking its
  correlation with treatment and outcome), imposes the empirical null on the
  outcome as usual, and builds one synthetic instrument per `rho` as
  `z_r = r * x_std + sqrt(1 - r^2) * L_perp_std`, whose in-sample Pearson
  correlation with the treatment is *exactly* `r` on every rep. The same
  recombined outcome and permuted `L` are shared across `rho` values (common
  random numbers), so the returned family of power curves isolates what
  instrument strength alone does to power. `print()` shows an MDE row and an
  implied first-stage F per `rho`; `plot()` draws one colour-coded curve per
  `rho`; `mde()` and `power_at()` return per-`rho` data frames.
  The default `z =` path is unchanged (`z` is now optional in the signature,
  but exactly one of `z` / `L` must be supplied). `plan_mde()` declines
  hypothetical-instrument objects with a clear error.

# simpower 0.1.4

* **RDD diagnostics now use consistent units.** When `id` is supplied,
  `power_rdd()`'s sample-composition check compared a *cluster* count (the
  treated-units target) against an *observation* count (the simulation
  average), which could read as a large mismatch (e.g. 50 vs 317 on the Lee
  senate data) even though both numbers were correct. The simulated treated
  count is now the number of treated clusters, the metric labels say
  explicitly what is counted and that the simulation average is measured
  inside the bandwidth window, and the bandwidth table's `n_in_bw` column now
  always counts observations in the window (previously clusters when `id` was
  given). Documented in `?power_rdd`.

# simpower 0.1.3

* **Fixed: `plot()` errored for event-study and RDD objects in a fresh
  session.** With the default `type = "both"`, the two panels were combined
  with `p1 + p2 + patchwork::plot_annotation(...)`. The `+` method that lets
  two ggplots be combined is only registered when the patchwork namespace is
  loaded, and the first `+` was evaluated *before* `patchwork::` triggered
  that load -- so `plot()` failed with "Can't add `.gg_horizon(x)` to a
  <ggplot> object" unless the user had loaded patchwork themselves. The
  panels are now combined with `patchwork::wrap_plots()`.
* **`mde_peak` added to `globalVariables()`**, silencing the R CMD check NOTE
  introduced with the per-horizon plot in 0.1.2.
* **Tests skip gracefully when `AER` is absent** (it is only in Suggests).
* **Maintainer email filled in** (was a placeholder, which blocks CRAN
  submission).

# simpower 0.1.2

* **Event-study horizon cap raised from 10 to 20.** `power_event()` now models
  up to 20 post-treatment horizons by default (still overridable with
  `horizons`).
* **Per-horizon MDE now takes the effect `shape` into account.** The per-horizon
  table gains an `mde_peak` column (`= mde / weight`) that expresses each
  horizon's MDE in units of the peak effect, and `plot(pe, type = "horizon")`
  uses it. For a `"constant"` shape (all weights 1) `mde_peak == mde`, so the
  graph is unchanged.
* **`plan_mde()` plot** now places x-axis ticks exactly on the sampled
  `n_units` / `n_treated`, so points line up with their labels.
* **References updated** (Doucette 2025, APSR; Doucette 2026, arXiv) and package
  metadata (author, repository URL).

# simpower 0.1.1

Correctness and robustness fixes from a code review.

* **Controls with missing values no longer error.** `.controls_matrix()` now
  builds its model frame with `na.pass`, so rows with `NA` controls survive to
  the common `complete.cases()` step instead of being silently dropped by
  `model.matrix()` (which caused a row-length-mismatch error in every design).
* **Non-numeric outcomes / regressors are rejected, not silently coerced.**
  Outcome, treatment, instrument and running-variable columns are pulled with a
  new `.numcol()` helper that errors on a factor/character column rather than
  turning it into integer level codes.
* **Event study no longer contaminates the baseline.** Observations at event
  times outside the estimation window (leads beyond `leads`, horizons beyond the
  cap) are now dropped with a message instead of being folded into the
  never-treated comparison group.
* **`crit_method` argument (default `"empirical"`).** The test's critical value
  is read off the simulated null distribution of the t-statistics, giving exact
  size under clustered / small-sample SEs (the old normal approximation was
  visibly under-sized). `crit_method = "normal"` restores the previous behaviour.
* **`use_fixest` now does something.** `power_twfe()` reports a real-data point
  estimate in `extras$estimate` (via `fixest::feols()` when available, else the
  internal within-estimator), matching the long-standing documentation.
* **Recombination uses a derangement**, so no unit is ever paired with its own
  outcome series (a plain permutation left ~1 fixed point per rep).
* **Degrees-of-freedom guard.** Rank-deficient / too-small subsamples now return
  an invalid fit (dropped) instead of a silent SE of 0, and a warning is issued
  when a meaningful share of reps is dropped.
* **`plan_mde()` robustness.** It re-runs against the data stored on the fitted
  object (no longer breaks when the original `data` binding is gone) and only
  builds its ggplot when `plot = TRUE`.
* **Performance.** Event-study fits demean the panel once (not twice) per rep;
  `.demean_twoway()` short-circuits balanced panels; `.subsample_units()` returns
  early on the full-sample path; `plan_mde()` skips the unused RDD bandwidth
  sweep; `sim_panel()`'s AR(1) errors are vectorised.
* **Packaging.** `AER` and `sandwich` moved to `Suggests` (diagnostic-only, now
  guarded by `requireNamespace()`); removed the spurious `LazyData: true`.

# simpower 0.1.0

* Initial release.
* Standardized MDE: every result also reports the MDE in within-SD-of-outcome
  units; `mde(x, standardized = TRUE)` and a "(SD units)" line in `print()`.
* `plan_mde()`: sweep `n_units` / `n_treated` and plot MDE vs sample size for
  ex-ante design planning.
* Query & tidy helpers: `power_at(x, effect)`, `mde(x, target = ...)` for
  arbitrary power levels, and `as.data.frame(x, what = ...)`.
* `power_rdd(kernel = "uniform")` in addition to the default triangular kernel.
* A warning is now issued when a user-supplied `emax` is too small for the power
  curve to reach a requested target (MDE would otherwise be a silent `NA`).
* Simulation-based power / MDE for four designs: `power_twfe()`,
  `power_event()`, `power_iv()`, `power_rdd()`.
* Common controls, `n_units`, and `n_treated` design-planning arguments.
* Event study: constant / increasing / decreasing effect shapes and a
  per-horizon MDE graph.
* RDD: adjustable bandwidth and an MDE-versus-bandwidth graph.
* `print()`, `summary()`, and `plot()` methods; ggplot2 figures.
