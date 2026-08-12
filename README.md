# simpower

**Simulation-based statistical power for causal designs — computed from your own data.**

`simpower` estimates statistical power and the **minimum detectable effect (MDE)**
for four causal-inference designs, using the Monte-Carlo *recombination* method
of Doucette (2025; 2026):

| Function | Design | Extra figure |
|---|---|---|
| `power_twfe()` | Two-way (unit + time) fixed effects | — |
| `power_event()` | Event study | treatment-by-post-period MDE graph |
| `power_iv()` | Instrumental variables (2SLS) | power curves by instrument strength (optional, via `L`/`rho`) |
| `power_rdd()` | Sharp regression discontinuity | MDE-vs-bandwidth graph |

## The approach

For each design the outcome for each unit is repeatedly randomly recombined with another units treatment so that any true
relationship to the treatment is broken — this manufactures the empirical null
while preserving the real features of the data that drive power
(e.g. autocorrelation, unbalanced panels, etc). Each design's estimator and standard errors are
computed on each recombined sample and stored. A growing true effect is then injected:
adding an effect `τ` shifts the estimated coefficient by exactly `τ` and leaves
the standard error unchanged. The whole power curve is then a cheap
post-processing step over one set of simulations, and the MDE is the smallest
effect whose power reaches the target (80% by default).

## What each design estimates

Every design fits a linear-in-outcome estimator on the recombined sample; the
**tested coefficient** is the one a growing effect is injected into.

**`power_twfe()` — two-way fixed effects.** Within regression, unit-clustered SE:

$$Y_{it} = \beta\,X_{it} + \gamma' Z_{it} + \alpha_i + \delta_t + \varepsilon_{it}$$

The coefficient on the treatment $X_{it}$, $\beta$, is tested. $\alpha_i$ and
$\delta_t$ are unit and time fixed effects, $Z_{it}$ optional controls.

**`power_event()` — event study.** Two regressions on the same panel. The
**overall** curve uses a single shaped-exposure term $P_{it}$ — the
post-treatment event-time dummies collapsed with the assumed `shape`
(`"constant"`/`"increasing"`/`"decreasing"`), normalised so the peak horizon has
weight 1 (so the MDE is in units of the *peak* effect):

$$Y_{it} = \beta\,P_{it} + \gamma' Z_{it} + \alpha_i + \delta_t + \varepsilon_{it}$$

The **treatment-by-post-period** graph comes from the dynamic regression on the event-time
dummies $D_{it}^{k}$ (relative time $k = t - \text{adoption}$; horizon $-1$
omitted as the baseline, leads kept as placebos):

$$Y_{it} = \sum_{k} \theta_k\,D_{it}^{k} + \gamma' Z_{it} + \alpha_i + \delta_t + \varepsilon_{it}$$

each $\theta_k$ giving the effect detectable at that time period; the treatment-by-post-period
graph reports these in units of the peak effect ($\theta_k$'s MDE divided by the
shape weight $w_k$), so it reflects the assumed `shape`. Both use unit-clustered
SEs.

**`power_iv()` — instrumental variables (2SLS).** Excluded instrument $Z_i$,
exogenous controls $W_i$:

$$Y_i = \beta\,D_i + \gamma' W_i + \varepsilon_i \qquad\text{(structural)}$$

$$D_i = \pi\,Z_i + \gamma' W_i + v_i \qquad\text{(first stage)}$$

The coefficient on the endogenous $D_i$, $\beta$, is tested; the first-stage
$F$ is reported. HC1 SE, or clustered when an `id` is supplied.

**`power_rdd()` — sharp regression discontinuity.** Local-linear fit on
observations within the bandwidth, kernel-weighted, with
$T_i = \mathbf{1}\{R_i \ge c\}$:

$$Y_i = \tau\,T_i + \beta_0 + \beta_1 (R_i - c) + \beta_2\,T_i (R_i - c) + \gamma' W_i + \varepsilon_i$$

The jump at the cutoff, $\tau$, is tested. Triangular kernel by default (or
uniform); HC1 SE, or clustered when an `id` is supplied.

## Installation

```r
# install.packages("remotes")
remotes::install_github("faranordk/simpower")
```

Required: **ggplot2**, **patchwork**, **AER**, **sandwich**, **lmtest**.
Optional accelerators used automatically if present: **fixest** (TWFE / event
study point estimates) and **rdrobust** (MSE-optimal RDD bandwidth selection).

## Shared arguments

Across all four designs you can set:

* `controls` — a character vector of control-variable columns to partial out;
* `n_units` — the number of units to draw for each simulated sample;
* `n_treated` — how many of those units are treated;
* `reps` (default `1000`) — the number of recombination simulations. More reps
  give a smoother power curve and a more precise MDE; they do **not** change the
  answer in expectation (no bias), only the Monte-Carlo noise around it. Start
  low (e.g. `200`) while exploring, then raise to `2000`–`5000` for a final run.
* `alpha` (default `0.05`) — the significance level of the test. Power is the
  share of simulated estimates that reach significance at this level, so a
  smaller `alpha` (a stricter test) means a larger MDE.
* `alternative` (default `"greater"`) — the direction of the test.
  `"greater"` / `"less"` are one-sided (the paper's convention: is the effect
  positive / negative?); `"two.sided"` tests either direction and, being
  stricter, yields a somewhat larger MDE.
* `power_targets` (default `c(0.8, 0.9)`) — the power level(s) at which the MDE
  is reported. `c(0.8, 0.9)` reports the smallest effect detectable 80% and 90%
  of the time; pass any values in `(0, 1)`, e.g. `power_targets = 0.8`.
* `seed` (default `5000`) — the random seed, so a run is exactly reproducible.

(A further option, `crit_method`, defaults to `"empirical"` — the critical
value is read off the simulated null t-statistics, giving the test exactly the
requested `alpha` even under clustered / small-sample SEs; `"normal"` uses the
Gaussian quantile instead.)

Every function returns a `simpower` object with `print()`, `summary()`, and
`plot()` methods and an `mde()` extractor.

Two further options are shared but rarely touched: `emax` and `step` set the
largest effect and the step size of the grid the power curve is evaluated on.
Left `NULL` (the default) they are chosen automatically from the null SE scale;
widen `emax` if a warning says the curve never reaches your target power.

## Design-specific arguments

Each function also takes arguments unique to its design — first the columns that
define the estimate, then design-specific controls.

### `power_twfe()`

* `y`, `x`, `id`, `time` — column names for the outcome, the treatment whose
  coefficient is tested, the unit id, and the time variable.
* `use_fixest` (default `TRUE`) — if the **fixest** package is installed, use
  `fixest::feols()` for the real-data point estimate reported in the printout;
  otherwise an internal within-estimator (which matches `feols`) is used. The
  simulation itself always uses the internal estimator, so this only affects the
  reported `real-data estimate` line, not the power curve.

### `power_event()`

* `y`, `id`, `time` — outcome, unit id, and time columns.
* `treat` — an **absorbing** 0/1 treatment indicator (once it turns on it stays
  on). Optional if `cohort` is supplied.
* `cohort` — an alternative to `treat`: a column giving each unit's adoption
  time on the same scale as `time`; never-treated units take `NA`/`Inf`.
  Relative event time is `time − adoption`.
* `shape` (default `"constant"`) — the assumed shape of the true effect path
  across post-treatment horizons `k = 0, 1, …, H`. The overall curve regresses
  the outcome on a single *shaped-exposure* term whose weight at horizon `k` is
  `w_k`, normalised so the peak horizon has weight 1 — so the tested coefficient
  is the **peak-horizon effect** and the effect at horizon `k` is `β · w_k`:
  * `"constant"` — `w_k = 1`: the same effect at every horizon (a level shift).
  * `"increasing"` — `w_k = (k+1)/(H+1)`: the effect ramps up linearly, smallest
    just after treatment and peaking at the last horizon `H`.
  * `"decreasing"` — `w_k = (H+1−k)/(H+1)`: the effect is largest immediately
    after treatment and fades linearly.

  Because the shape decides whether identifying signal sits at early
  well-measured post-treatment time periods or later noisier ones, it moves the overall MDE, and
  the treatment-by-post-period graph now reflects it too (see below).

  The path is always **linear**, so with `"increasing"` the effect grows by a
  *constant* increment of `1/(H+1)` of the peak each period — `β/(H+1)` per
  horizon, not an accelerating curve. The **initial** effect (at `k = 0`, the
  period treatment turns on) is `w_0 = 1/(H+1)` of the peak, i.e. `β/(H+1)`; its
  size is set entirely by `H+1` = the number of post-treatment horizons (the
  `horizons` argument). For example with `H = 20` the effect starts at
  `1/21 ≈ 0.048` of the peak and adds another `1/21` each period up to `β` at
  horizon 20; setting `horizons = 4` makes `H+1 = 5`, so the effect starts at
  `1/5 = 0.2` of the peak and climbs in steps of `0.2β`. `"decreasing"` is the
  mirror image: it starts at the peak `β` at horizon 0 and falls by `β/(H+1)`
  each period.
* `horizons` — **how many periods after treatment to model**, i.e. the maximum
  post-treatment horizon `H` (default: as many as the data support, capped at
  20). This is the post-treatment counterpart of `leads`: `leads` sets how many
  periods *before* treatment enter as placebos, `horizons` sets how many periods
  *after* treatment are estimated. It also fixes the shape denominator `H+1`
  above, so shrinking it makes an `"increasing"`/`"decreasing"` path steeper.
* `leads` — the number of pre-treatment leads kept as placebo terms (default up
  to 5). Horizon `−1` is always the omitted baseline; event times outside the
  `[−leads, horizons]` window are dropped rather than folded into the
  never-treated comparison.

The **per-horizon graph** takes the shape into account. Each horizon's MDE is
divided by its weight `w_k` and reported in units of the *peak* effect
(`mde_peak = mde / w_k`), matching the overall MDE. Under `"increasing"`, early
horizons carry little assumed effect (small `w_k`), so a larger peak effect is
needed before they become individually detectable — the bars are tallest at the
early horizons and shrink toward the peak. For `"constant"` all `w_k = 1`, so
`mde_peak` equals the raw per-horizon MDE and the graph is unchanged.

### `power_iv()`

* `y`, `d`, `z` — the outcome, the **endogenous** regressor whose coefficient is
  tested, and the **excluded instrument**.
* `id` (optional) — a cluster id; when supplied, SEs are clustered by `id`
  instead of the default heteroskedasticity-robust (HC1).
* `time` (optional) — a time column, if the data are panel-structured.
* `L`, `rho` (optional) — switch on **hypothetical-instrument mode** (supply
  either `z` or `L`, not both). Use this *before you have an instrument*: it
  answers "if I could find an instrument with a first-stage correlation of r,
  how powered would this design be?". `L` is any variable in your data whose
  distributional character the imagined instrument should inherit (its scale
  and shape are used as the noise carrier; its actual correlations are not).
  On every simulation rep the `L` series is randomly permuted — breaking its
  correlation with both the treatment and the outcome — and one synthetic
  instrument is built per value of `rho` as

  $$z_r = r\,\tilde X + \sqrt{1-r^2}\,\tilde L_\perp,$$

  where $\tilde X$ is the standardised treatment and $\tilde L_\perp$ the
  permuted `L` residualised on the treatment and standardised. By construction
  the in-sample Pearson correlation between $z_r$ and the treatment is
  **exactly** `r` on every rep, and the same recombined outcome is shared
  across the `rho` values — so the returned family of power curves isolates
  what instrument strength alone does to power. `rho` defaults to
  `c(0.1, 0.3, 0.5)`; each value must be strictly between 0 and 1 (the sign
  is irrelevant: 2SLS is invariant to the instrument's sign). The printout
  reports one MDE row per `rho` plus the implied first-stage F
  ($F = r^2 (n-2) / (1-r^2)$, deterministic because the correlation is exact),
  and `plot()` draws one colour-coded curve per `rho`. This mode is
  cross-sectional (omit `time`; `id` may still be supplied for clustered SEs),
  and `mde()` / `power_at()` return per-`rho` data frames.

### `power_rdd()`

* `y`, `run` — the outcome and the running (forcing) variable.
* `cutoff` (default `0`) — the threshold of `run` at which treatment switches on
  (`treatment = run >= cutoff`).
* `bandwidth` — the half-width of the estimation window around the cutoff. When
  `NULL`, an MSE-optimal bandwidth from **rdrobust** is used if available,
  otherwise a rule-of-thumb value (with a message).
* `kernel` (default `"triangular"`) — the weighting kernel for the local-linear
  fit; `"uniform"` is equivalent to OLS within the window.
* `bw_grid` — the vector of bandwidths for the MDE-vs-bandwidth graph (default:
  a grid spanning roughly 0.5×–1.5× the main bandwidth).
* `bw_reps` (default `300`) — simulations per bandwidth in that sweep, kept
  smaller than `reps` to bound run time.
* `id` (optional) — a cluster id; SEs are clustered by `id` when supplied.
* `use_rdrobust` (default `TRUE`) — use **rdrobust** for MSE-optimal bandwidth
  selection when `bandwidth = NULL`.

## Examples

Each block uses a built-in data generator (`sim_*`) so it runs as-is; swap in
your own data frame and column names to use it for real. Console output is
illustrative — the exact MDE depends on your data.

### 1. TWFE

```r
library(simpower)

df <- sim_panel(n_units = 100, n_time = 10)
pw <- power_twfe(df, y = "y", x = "x", id = "id", time = "t",
                 controls = NULL, n_units = 100, n_treated = 50,
                 reps = 500)
pw
```

```
<simpower> Two-way fixed effects
  test: one-sided (greater), alpha = 0.05, reps = 500
  null: mean beta = 0.0003, mean se = 0.0611
  MDE:  0.8 power = 0.152   0.9 power = 0.179
  MDE (SD units):  0.8 power = 0.28   0.9 power = 0.33   [within-SD of Y = 0.54]
  real-data estimate: b = 0.306 (se = 0.058, within)
  use summary() for the sample-composition check and plot() for figures.
```

```r
summary(pw)   # + sample-composition check (target vs simulated draw)
plot(pw)      # the power curve

## variants -----------------------------------------------------------------
# two-sided test
power_twfe(df, "y", "x", "id", "t", alternative = "two.sided", reps = 500)

# plan a smaller design: 40 units, 20 treated
power_twfe(df, "y", "x", "id", "t", n_units = 40, n_treated = 20, reps = 500)

# with your own data: partial out controls
# power_twfe(mydata, "y", "x", "id", "t", controls = c("z1", "z2"), reps = 500)
```

### 2. Event study

```r
ev <- sim_event(n_units = 120, n_time = 12, adopt = 7)
pe <- power_event(ev, y = "y", id = "id", time = "t", treat = "d",
                  shape = "increasing", reps = 500)   # or "constant" / "decreasing"
pe
```

```
<simpower> Event study
  test: one-sided (greater), alpha = 0.05, reps = 500
  null: mean beta = 0.0011, mean se = 0.0740
  MDE:  0.8 power = 0.184   0.9 power = 0.217   (MDE in units of the peak effect)
  effect shape: increasing
  use summary() for the per-horizon MDE table and plot() for figures.
```

```r
plot(pe, type = "overall")    # overall power curve
plot(pe, type = "horizon")    # per-horizon MDE (which horizons are well powered)
plot(pe)                      # both, side by side (needs current patchwork/ggplot2)
summary(pe)                   # prints the per-horizon MDE table

## cap the horizons/leads and try a constant shape
power_event(ev, y = "y", id = "id", time = "t", treat = "d",
            shape = "constant", horizons = 6, leads = 4, reps = 500)

# with your own data: supply a cohort/adoption-time column instead of a 0/1 indicator
# power_event(mydata, "y", "id", "t", cohort = "adopt_year", shape = "constant")
```

### 3. Instrumental variables (2SLS)

```r
iv <- sim_iv(n = 1000, strength = 0.4)
pv <- power_iv(iv, y = "y", d = "d", z = "z", reps = 500)
pv
```

```
<simpower> Instrumental variables (2SLS)
  test: one-sided (greater), alpha = 0.05, reps = 500
  null: mean beta = 0.0009, mean se = 0.1032
  MDE:  0.8 power = 0.257   0.9 power = 0.302
  first-stage F: 43.9
```

```r
plot(pv)                      # the power curve

## weak vs strong instrument: watch the MDE grow as strength falls
power_iv(sim_iv(n = 1000, strength = 0.15), y = "y", d = "d", z = "z", reps = 500)

# with your own data: add exogenous controls and cluster SEs by id
# power_iv(mydata, "y", "d", "z", controls = c("w1", "w2"), id = "cluster")
```

#### Hypothetical instrument: what strength would I need?

No instrument yet? Supply `L` instead of `z` and a set of assumed first-stage
correlations, and see how the power curve shifts with instrument strength:

```r
iv <- sim_iv(n = 1000, strength = 0.4)
iv$cand <- rgamma(nrow(iv), 2)     # any series; only its scale/shape are used
ph <- power_iv(iv, y = "y", d = "d", L = "cand",
               rho = c(0.1, 0.3, 0.5), reps = 500)
ph
```

```
<simpower> Instrumental variables (2SLS)
  test: one-sided (greater), alpha = 0.05, reps = 500
  hypothetical instrument: L = 'cand', curves by cor(z, x)
  MDE by instrument strength:
 rho   80%   90% med_se implied_F 80% (SD) 90% (SD)
 0.1 0.907 1.102 0.3902      10.1    0.485    0.589
 0.3 0.322 0.371 0.1266      98.7    0.172    0.198
 0.5 0.195 0.234 0.0758     332.7    0.104    0.125
```

```r
plot(ph)                      # one colour-coded power curve per rho
mde(ph, target = 0.7)         # per-rho MDE at any power level (data frame)
power_at(ph, c(0.25, 0.5))    # power at given effects, one column per rho
```

Reading it: with n = 1000, an instrument correlating only 0.1 with the
treatment (implied first-stage F ≈ 10) needs a true effect of ~0.91 for 80%
power, while r = 0.5 detects ~0.20 — a direct answer to "how strong an
instrument do I need to detect the effect I expect?". Each rep permutes `L`
(so it carries no real correlation with anything) and rebuilds the synthetic
instrument at exactly the requested correlation; see the argument description
above for the construction.

### 4. Regression discontinuity

```r
rd <- power_rdd(sim_rdd(n = 2000), y = "y", run = "x", cutoff = 0,
                bandwidth = 0.5, reps = 500)
rd
```

```
<simpower> Regression discontinuity
  test: one-sided (greater), alpha = 0.05, reps = 500
  null: mean beta = -0.0021, mean se = 0.1187
  MDE:  0.8 power = 0.295   0.9 power = 0.347
  bandwidth: 0.5  (user-supplied)
```

```r
plot(rd, type = "overall")    # power curve
plot(rd, type = "bandwidth")  # MDE vs bandwidth

# let rdrobust pick an MSE-optimal bandwidth (bandwidth = NULL)
power_rdd(sim_rdd(n = 2000), y = "y", run = "x", cutoff = 0, reps = 500)

# uniform kernel (= OLS within the window); custom bandwidth sweep
power_rdd(sim_rdd(n = 2000), y = "y", run = "x", cutoff = 0, bandwidth = 0.5,
          kernel = "uniform", bw_grid = c(0.25, 0.5, 0.75, 1.0), reps = 500)
```

> **Plotting note.** For event studies and RDD, `plot(x)` defaults to
> `type = "both"`, which stitches the two panels together with **patchwork**
> (this works in a fresh session as of 0.1.3 — no `library(patchwork)`
> needed). Draw one panel at a time with `type = "overall"` /
> `"bandwidth"` (`"horizon"` for event studies).

## Reading the result

```r
mde(pw)                          # MDE at the fitted power targets
mde(pw, target = 0.7)            # MDE at an arbitrary power level
mde(pw, standardized = TRUE)     # MDE in within-SD-of-outcome units
power_at(pw, c(0.1, 0.2))        # power at given true effect sizes
as.data.frame(pw)                # the power curve as a data frame
```

Every `print()` also shows the MDE in standard-deviation units (raw MDE divided
by the within-unit SD of the outcome), so designs with different outcomes are
directly comparable.

## Planning helpers

```r
## how many units do I need? -- sweep sample size and plot MDE vs N
plan_mde(pw, n_units = c(50, 100, 200, 400))
plan_mde(pw, n_treated = c(10, 25, 50))       # or sweep the treated count
```

## Check how power depends on design choices

* **More units / more treated units** shrink the MDE (more clusters, more
  identifying variation). Vary `n_units` / `n_treated` to plan sample size.
* **Event-study shape**: an `"increasing"` effect concentrates identifying
  signal at later, noisier horizons, so its overall MDE (expressed in units of
  the *peak* effect) differs from a `"constant"` one; the per-horizon graph
  shows exactly which horizons are well powered.
* **RDD bandwidth**: a wider bandwidth uses more data (smaller MDE) at the cost
  of bias; the MDE-vs-bandwidth graph makes that trade-off explicit.
* **IV instrument strength**: in hypothetical-instrument mode (`L` + `rho`),
  the fan of power curves shows directly how much power a stronger first stage
  buys — useful for deciding whether a candidate instrument is worth pursuing
  before committing to a design.

## Notes on standard errors

* TWFE / event study: unit-clustered SEs with the Stata/`reghdfe`
  small-sample factor.
* IV: heteroskedasticity-robust (HC1), or clustered when an `id` is supplied.
* RDD: local-linear with a triangular kernel and HC1 (or clustered) SEs.

## References

Doucette, M. (2025). What Can We Learn About the Effects of Democracy Using
Cross-National Data? *American Political Science Review*, 119(3).

Doucette, M. (2026). simpower: Simulation-Based Statistical Power for Panel and
Causal Designs. *arXiv*.
