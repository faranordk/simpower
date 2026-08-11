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
