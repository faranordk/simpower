# Regression tests for the review fixes.

test_that("controls containing NA no longer crash (na.pass alignment)", {
  set.seed(101)
  df <- sim_panel(50, 8)
  df$z1 <- rnorm(nrow(df))
  df$z1[c(3, 17, 50)] <- NA                 # missing controls
  pw <- power_twfe(df, "y", "x", "id", "t", controls = "z1", reps = 80)
  expect_s3_class(pw, "simpower")
  expect_true(is.finite(pw$mde["80%"]))
})

test_that("factor / character outcome or treatment is rejected, not coerced", {
  set.seed(102)
  df <- sim_panel(40, 6)
  df$xf <- factor(ifelse(df$x > 0, "hi", "lo"))
  expect_error(power_twfe(df, "y", "xf", "id", "t", reps = 20),
               "must be numeric")
})

test_that("event study drops out-of-window horizons instead of baselining them", {
  set.seed(103)
  de <- sim_event(n_units = 60, n_time = 30, adopt = 4)  # true max horizon = 26
  expect_message(
    pw <- power_event(de, "y", "id", "t", treat = "d", reps = 60),
    "outside the window"
  )
  expect_equal(pw$extras$H, 20L)
  expect_true(all(pw$extras$Kpost <= 20))
})

test_that("empirical critical value gives ~nominal size; normal option still works", {
  set.seed(104)
  df <- sim_panel(60, 8)
  emp <- power_twfe(df, "y", "x", "id", "t", reps = 1500, crit_method = "empirical")
  nrm <- power_twfe(df, "y", "x", "id", "t", reps = 1500, crit_method = "normal")
  # size = power at effect 0 should sit right at alpha under the empirical calibration
  expect_lt(abs(emp$power$power[1] - 0.05), 0.012)
  # the two methods differ (empirical corrects the under-sizing)
  expect_false(isTRUE(all.equal(emp$power$power[1], nrm$power$power[1])))
})

test_that("power_from_null empirical calibration is exact on a synthetic null", {
  set.seed(105)
  bs <- rnorm(4000); ses <- rep(1, 4000)
  g  <- seq(0, 5, by = 0.05)
  p  <- power_from_null(bs, ses, g, alpha = 0.05, alternative = "greater",
                        crit_method = "empirical")
  expect_lt(abs(p[1] - 0.05), 0.005)          # size hit essentially exactly
  expect_true(all(diff(p) >= -1e-9))          # still monotone
})

test_that("recombination uses a derangement (no unit keeps its own series)", {
  set.seed(106)
  fp <- replicate(500, { d <- simpower:::.derangement(60); sum(d == seq_len(60)) })
  expect_equal(sum(fp), 0)                     # never a fixed point
  expect_true(all(vapply(1:50, function(i) {
    d <- simpower:::.derangement(40); setequal(d, 1:40)   # still a permutation
  }, logical(1))))
})

test_that(".ols_vcov signals an invalid fit (NULL) when dof are exhausted", {
  set.seed(107)
  X <- cbind(1, rnorm(6)); y <- rnorm(6)
  expect_null(simpower:::.ols_vcov(X, y, dfK = 8))    # dfK > n -> NULL, not SE 0
  fit <- simpower:::.ols_vcov(X, y, dfK = 2)          # healthy case still works
  expect_true(is.finite(fit$V[1, 1]))
})

test_that("combined event fit is numerically identical to the two separate fits", {
  set.seed(108)
  n <- 300
  uid <- rep(1:30, each = 10); ti <- rep(1:10, times = 30)
  y <- rnorm(n)
  P <- runif(n)
  D <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  comb <- simpower:::.fit_event_combined(y, P, D, uid, ti)
  ov   <- simpower:::.fit_twfe(y, P, uid, ti)
  dy   <- simpower:::.fit_event_dynamic(y, D, uid, ti)
  expect_equal(unname(comb$ov["b"]),  unname(ov["b"]),  tolerance = 1e-9)
  expect_equal(unname(comb$ov["se"]), unname(ov["se"]), tolerance = 1e-9)
  expect_equal(comb$b,  dy$b,  tolerance = 1e-9)
  expect_equal(comb$se, dy$se, tolerance = 1e-9)
})

test_that("TWFE now reports a real-data point estimate", {
  set.seed(109)
  df <- sim_panel(60, 8, beta = 0.4)
  pw <- power_twfe(df, "y", "x", "id", "t", reps = 60)
  expect_true(all(is.finite(pw$extras$estimate)))
  expect_true(!is.null(pw$extras$estimator))
})

test_that("plan_mde survives loss of the original data binding and plot = FALSE", {
  set.seed(110)
  local({
    d2 <- sim_panel(80, 8)
    pw <- power_twfe(d2, "y", "x", "id", "t", reps = 60)
    rm(d2)                                     # original binding gone
    pl <- plan_mde(pw, n_units = c(30, 60), reps = 50, plot = FALSE)
    expect_s3_class(pl, "simpower_plan")
    expect_equal(nrow(pl$table), 2)
    expect_null(pl$plot)                       # no ggplot built when plot = FALSE
  })
})

test_that("subsample early-return keeps the full sample when nothing requested", {
  uid <- rep(1:10, each = 4)
  ss <- simpower:::.subsample_units(uid, treated_ids = 1:5,
                                    n_units = NULL, n_treated = NULL)
  expect_identical(ss$rows, seq_along(uid))    # every row, original order
  expect_identical(ss$uid, uid)                # original unit ids untouched
})

test_that("two.sided and less run end-to-end on a real design", {
  set.seed(111)
  df <- sim_panel(50, 8)
  p2 <- power_twfe(df, "y", "x", "id", "t", reps = 120, alternative = "two.sided")
  pl <- power_twfe(df, "y", "x", "id", "t", reps = 120, alternative = "less")
  expect_true(is.finite(p2$mde["80%"]))
  expect_true(is.finite(pl$mde["80%"]))
  expect_lt(abs(p2$power$power[1] - 0.05), 0.02)
})

test_that("plot() composes event/RDD panels without patchwork attached", {
  set.seed(112)
  de <- sim_event(n_units = 50, n_time = 8, adopt = 5)
  pe <- power_event(de, "y", "id", "t", treat = "d", reps = 40)
  expect_s3_class(plot(pe), "patchwork")            # type = "both" (default)
  expect_s3_class(plot(pe, type = "overall"), "ggplot")
  dr <- sim_rdd(n = 800, jump = 0)
  pr <- power_rdd(dr, "y", "x", cutoff = 0, bandwidth = 0.5,
                  reps = 40, bw_reps = 30)
  expect_s3_class(plot(pr), "patchwork")
})

test_that("RDD diagnostics use consistent units (clusters vs observations)", {
  set.seed(113)
  dr <- sim_rdd(n = 600, jump = 0)
  dr$cl <- sample(1:30, nrow(dr), replace = TRUE)   # 30 clusters
  pr <- power_rdd(dr, "y", "x", cutoff = 0, bandwidth = 0.5, id = "cl",
                  reps = 40, bw_reps = 30)
  sc <- pr$sample_check
  ## with `id`, both rows count clusters -> sim averages bounded by 30
  expect_lte(sc$sim_avg[1], 30)
  expect_lte(sc$sim_avg[2], 30)
  expect_match(sc$metric[2], "treated units")
  ## n_in_bw counts observations, not clusters
  bt <- pr$extras$bandwidth
  expect_true(all(bt$n_in_bw > 30))
})

test_that("hypothetical-instrument mode: exact correlation and monotone power", {
  set.seed(114)
  ## .synth_instrument hits the requested correlation exactly
  x <- rnorm(300); Lp <- sample(rexp(300))
  for (r in c(0.1, 0.3, 0.5)) {
    z <- simpower:::.synth_instrument(x, Lp, r)
    expect_equal(cor(z, x), r, tolerance = 1e-12)
  }
  ## end-to-end: stronger instrument -> smaller MDE, higher power everywhere
  df <- sim_iv(n = 500, strength = 0.5)
  df$cand <- rgamma(nrow(df), 2)          # skewed candidate series
  ph <- power_iv(df, "y", "d", L = "cand", rho = c(0.1, 0.3, 0.5), reps = 300)
  expect_s3_class(ph, "simpower")
  expect_true(all(c("effect", "power", "rho") %in% names(ph$power)))
  m <- mde(ph)
  expect_true(is.data.frame(m) && nrow(m) == 3)
  expect_true(all(diff(m[["80%"]]) < 0))  # MDE falls as rho rises
  pa <- power_at(ph, c(0.5, 1))
  expect_true(all(pa$r0.5 >= pa$r0.1))    # stronger instrument dominates
  ## methods run
  expect_s3_class(plot(ph), "ggplot")
  expect_output(print(ph), "instrument strength")
  ## implied F: r^2 (n-2) / (1-r^2)
  expect_equal(unname(ph$extras$first_stage_F["r0.3"]),
               0.09 * (500 - 2) / 0.91, tolerance = 1e-6)
})

test_that("hypothetical-instrument mode argument validation", {
  df <- sim_iv(n = 200)
  df$cand <- rnorm(200)
  expect_error(power_iv(df, "y", "d"), "Supply `z`")
  expect_error(power_iv(df, "y", "d", z = "z", L = "cand"), "not both")
  expect_error(power_iv(df, "y", "d", L = "cand", rho = c(0, 0.5)),
               "strictly between 0 and 1")
  ph <- power_iv(df, "y", "d", L = "cand", rho = 0.4, reps = 60)
  expect_error(plan_mde(ph, n_units = c(50, 100)), "does not support")
  ## default path unchanged
  pw <- power_iv(df, "y", "d", z = "z", reps = 60)
  expect_true(is.numeric(pw$mde) && is.finite(pw$mde["80%"]))
})

## ---------------------------------------------------------------------------
## Control-preserving panel recombination (efficiency gains from controls)
## ---------------------------------------------------------------------------

test_that("recombine_panel_ctrl with no controls is byte-identical (same RNG stream)", {
  set.seed(107)
  df <- sim_panel(40, 6)
  uid <- as.integer(factor(df$id)); ti <- as.integer(factor(df$t))
  set.seed(1); a <- simpower:::.recombine_panel(df$y, uid, ti)
  set.seed(1); b <- simpower:::.recombine_panel_ctrl(df$y, uid, ti, W = NULL)
  expect_identical(a, b)
})

test_that("recombination preserves the controls' predicted component of y", {
  set.seed(108)
  df <- sim_panel(60, 8)
  df$w1 <- rnorm(nrow(df))
  df$y  <- df$y + 2 * df$w1                # w1 strongly predicts y within
  uid <- as.integer(factor(df$id)); ti <- as.integer(factor(df$t))
  W   <- matrix(df$w1, ncol = 1)
  set.seed(2); rec <- simpower:::.recombine_panel_ctrl(df$y, uid, ti, W)
  keep <- rec$keep
  ## the recombined outcome must still be predicted by the (own) controls...
  Md <- simpower:::.demean_twoway(cbind(W[keep, , drop = FALSE], rec$y[keep]),
                       as.integer(factor(uid[keep])), as.integer(factor(ti[keep])))
  g  <- coef(lm(Md[, 2] ~ Md[, 1] - 1))
  expect_equal(unname(g), 2, tolerance = 0.15)
  ## ...while the plain recombination destroys that relationship
  set.seed(2); rec0 <- simpower:::.recombine_panel(df$y, uid, ti)
  k0 <- rec0$keep
  Md0 <- simpower:::.demean_twoway(cbind(W[k0, , drop = FALSE], rec0$y[k0]),
                        as.integer(factor(uid[k0])), as.integer(factor(ti[k0])))
  g0 <- coef(lm(Md0[, 2] ~ Md0[, 1] - 1))
  expect_lt(abs(unname(g0)), 0.5)
})

test_that("predictive controls now shrink the simulated SE and the MDE", {
  set.seed(109)
  df <- sim_panel(80, 8)
  df$w1 <- rnorm(nrow(df))
  df$y  <- df$y + 2 * df$w1
  p0 <- power_twfe(df, "y", "x", "id", "t", reps = 400)
  p1 <- power_twfe(df, "y", "x", "id", "t", controls = "w1", reps = 400)
  ## controls explain a large share of within variance -> smaller SE and MDE
  expect_lt(p1$null$se_mean, 0.75 * p0$null$se_mean)
  expect_lt(p1$mde["80%"], 0.75 * p0$mde["80%"])
  ## null estimates stay centered at zero, spread shrinks with the SE
  expect_lt(abs(p1$null$b_mean), 0.05)
  expect_lt(sd(p1$null$bs), 0.75 * sd(p0$null$bs))
  ## size at effect 0 remains ~alpha under the empirical calibration
  expect_lt(abs(p1$power$power[1] - 0.05), 0.02)
})

test_that("pure-noise controls leave the answer essentially unchanged", {
  set.seed(110)
  df <- sim_panel(80, 8)
  df$junk <- rnorm(nrow(df))
  p0 <- power_twfe(df, "y", "x", "id", "t", reps = 300)
  p1 <- power_twfe(df, "y", "x", "id", "t", controls = "junk", reps = 300)
  expect_equal(p1$null$se_mean, p0$null$se_mean, tolerance = 0.05)
  expect_equal(unname(p1$mde["80%"]), unname(p0$mde["80%"]), tolerance = 0.08)
})

test_that("event study and panel IV use the control-preserving recombination", {
  set.seed(111)
  ev <- sim_event(n_units = 70, n_time = 10, adopt = 6)
  ev$w1 <- rnorm(nrow(ev))
  ev$y  <- ev$y + 2 * ev$w1
  e0 <- power_event(ev, "y", "id", "t", treat = "d", reps = 250)
  e1 <- power_event(ev, "y", "id", "t", treat = "d", controls = "w1", reps = 250)
  expect_lt(e1$null$se_mean, 0.75 * e0$null$se_mean)
})

test_that("panel IV uses an estimator-consistent (no-FE) control fit", {
  set.seed(112)
  n_units <- 60; n_time <- 6; n <- n_units * n_time
  id <- rep(1:n_units, each = n_time); t <- rep(1:n_time, n_units)
  w  <- rnorm(n)                                # control predicting y in levels
  z  <- rnorm(n)                                # instrument
  dd <- 0.7 * z + rnorm(n)                      # first stage
  y  <- 2 * w + rnorm(n)                        # null: no effect of dd
  df <- data.frame(y, d = dd, z, w, id, t)
  p0 <- power_iv(df, "y", "d", z = "z", id = "id", time = "t", reps = 300)
  p1 <- power_iv(df, "y", "d", z = "z", id = "id", time = "t",
                 controls = "w", reps = 300)
  ## strong control -> smaller SE; and never larger than without controls
  expect_lt(p1$null$se_mean, 0.75 * p0$null$se_mean)
  expect_lt(abs(p1$null$b_mean), 0.1)
})

## ---------------------------------------------------------------------------
## Upsampling: n_units / n_treated beyond the sample (cluster bootstrap)
## ---------------------------------------------------------------------------

test_that(".subsample_units upsamples with full copies + remainder, distinct instances", {
  uid <- rep(1:10, each = 4)
  set.seed(201)
  ss <- simpower:::.subsample_units(uid, treated_ids = 1:5, n_units = 25)
  expect_equal(length(unique(ss$uid)), 25)          # 25 distinct instances
  expect_equal(length(ss$rows), 25 * 4)             # each with its full series
  ## every original unit appears at least floor(25/10) = 2 times
  orig <- rep(uid[!duplicated(uid)], 0)             # noop; check via rows
  drawn <- uid[ss$rows[!duplicated(ss$uid[order(ss$rows)])]]
  expect_true(all(table(uid[ss$rows]) >= 2 * 4))
  ## within-sample request still returns plain ascending rows, original ids
  set.seed(202)
  ss2 <- simpower:::.subsample_units(uid, treated_ids = 1:5, n_units = 6)
  expect_true(!is.unsorted(ss2$rows))
  expect_equal(length(unique(ss2$uid)), 6)
  expect_true(all(ss2$uid %in% uid))
})

test_that("upsampling n_units shrinks the SE and the MDE (more data helps)", {
  set.seed(203)
  df <- sim_panel(50, 8)
  p1 <- power_twfe(df, "y", "x", "id", "t", reps = 300)
  p2 <- power_twfe(df, "y", "x", "id", "t", n_units = 100, reps = 300)
  ## doubling the clusters: SE should fall by roughly 1/sqrt(2)
  expect_lt(p2$null$se_mean, 0.8 * p1$null$se_mean)
  expect_lt(p2$mde["80%"], 0.85 * p1$mde["80%"])
  ## the sample-composition check hits the upsampled target
  expect_equal(p2$sample_check$sim_avg[1], 100, tolerance = 0.02)
})

test_that("n_treated beyond the treated pool keeps all controls and grows the total", {
  set.seed(204)
  df <- sim_panel(60, 6)                      # 30 treated, 30 control
  n_ctl <- 60 - length(unique(df$id[ave(df$x, df$id, FUN=function(v) diff(range(v)))>0]))
  p <- power_twfe(df, "y", "x", "id", "t", n_treated = 60, reps = 200)
  expect_equal(p$sample_check$sim_avg[2], 60, tolerance = 0.02)   # treated target met
  expect_gte(p$sample_check$sim_avg[1], 60 + n_ctl - 1)           # controls kept
})

test_that("upsampling works in the other designs (smoke)", {
  set.seed(205)
  iv <- sim_iv(n = 300)
  pv <- power_iv(iv, "y", "d", z = "z", n_units = 600, reps = 150)
  expect_true(is.finite(pv$null$se_mean))
  expect_equal(mean(pv$null$ses) < 1.5, TRUE)
  rd <- sim_rdd(n = 500)
  pr <- power_rdd(rd, "y", run = "x", cutoff = 0, bandwidth = 0.5,
                  n_units = 1000, reps = 150, bw_reps = 30)
  expect_true(is.finite(pr$null$se_mean))
  ev <- sim_event(n_units = 50, n_time = 8, adopt = 5)
  pe <- power_event(ev, "y", "id", "t", treat = "d", n_units = 100,
                    reps = 150)
  expect_true(is.finite(pe$null$se_mean))
})

## ---------------------------------------------------------------------------
## Heterogeneous effects (the `het` option)
## ---------------------------------------------------------------------------

test_that("het injection formulas are EXACT vs literally refitting (twfe)", {
  set.seed(301)
  df <- sim_panel(30, 6)
  uid <- as.integer(factor(df$id)); ti <- as.integer(factor(df$t))
  m <- sample(c(0, 2), 30, replace = TRUE)[uid]
  d <- m * df$x
  base <- simpower:::.fit_twfe(df$y, df$x, uid, ti, inject = d)
  for (e in c(0.3, 1.1, -0.7)) {
    ref <- simpower:::.fit_twfe(df$y + e * d, df$x, uid, ti)
    expect_equal(unname(base["b"] + e * base["lam"]), unname(ref["b"]), tolerance = 1e-10)
    v <- base["se"]^2 + e * base["v1"] + e^2 * base["v2"]
    expect_equal(unname(sqrt(v)), unname(ref["se"]), tolerance = 1e-10)
  }
})

test_that("het injection formulas are EXACT vs refitting (2SLS and RDD)", {
  set.seed(302)
  df <- sim_iv(n = 300)
  m <- sample(c(0.5, 1.5), 300, replace = TRUE)
  d <- m * df$d
  base <- simpower:::.fit_tsls(df$y, df$d, matrix(df$z, ncol = 1), inject = d)
  for (e in c(0.4, -0.9)) {
    ref <- simpower:::.fit_tsls(df$y + e * d, df$d, matrix(df$z, ncol = 1))
    expect_equal(unname(base["b"] + e * base["lam"]), unname(ref["b"]), tolerance = 1e-9)
    expect_equal(unname(sqrt(base["se"]^2 + e * base["v1"] + e^2 * base["v2"])),
                 unname(ref["se"]), tolerance = 1e-9)
  }
  rd <- sim_rdd(n = 800)
  rc <- rd$x; tr <- as.numeric(rc >= 0); w <- pmax(0, 1 - abs(rc) / 0.6)
  keep <- abs(rc) <= 0.6
  m2 <- sample(c(0, 2), sum(keep), replace = TRUE)
  d2 <- m2 * tr[keep]
  base <- simpower:::.fit_rdd(rd$y[keep], rc[keep], tr[keep], w[keep], inject = d2)
  for (e in c(0.5, 1.3)) {
    ref <- simpower:::.fit_rdd(rd$y[keep] + e * d2, rc[keep], tr[keep], w[keep])
    expect_equal(unname(base["b"] + e * base["lam"]), unname(ref["b"]), tolerance = 1e-9)
    expect_equal(unname(sqrt(base["se"]^2 + e * base["v1"] + e^2 * base["v2"])),
                 unname(ref["se"]), tolerance = 1e-9)
  }
})

test_that("het leaves the null draws bit-identical and size at alpha", {
  set.seed(303)
  df <- sim_panel(50, 8)
  p0 <- power_twfe(df, "y", "x", "id", "t", reps = 250)
  p1 <- power_twfe(df, "y", "x", "id", "t", reps = 250, het = c(0, 2))
  expect_identical(p0$null$bs, p1$null$bs)       # isolated het RNG stream
  expect_identical(p0$null$ses, p1$null$ses)
  ## at effect 0 the het and constant setups coincide -> same realized size
  expect_equal(p1$power$power[1], p0$power$power[1], tolerance = 1e-12)
  ## degenerate multiset = constant effects exactly (same explicit grid)
  p2 <- power_twfe(df, "y", "x", "id", "t", reps = 250, het = c(3, 3),
                   emax = max(p0$grid), step = p0$grid[2] - p0$grid[1])
  expect_equal(p0$power$power, p2$power$power, tolerance = 1e-12)
  ## heterogeneity costs power at the target: MDE weakly larger
  expect_gte(p1$mde["80%"], p0$mde["80%"])
  ## constant-effect comparison is stored and matches the plain run (up to the
  ## widened het grid's step size)
  step1 <- p1$grid[2] - p1$grid[1]
  expect_lt(abs(p1$extras$het$mde_constant["80%"] - p0$mde["80%"]), 2 * step1)
})

test_that("het spec validation and normalisation", {
  df <- sim_panel(30, 5)
  expect_error(power_twfe(df, "y", "x", "id", "t", reps = 30, het = c(-2, 0)),
               "positive mean")
  expect_error(power_twfe(df, "y", "x", "id", "t", reps = 30, het = "half"),
               "must be NULL")
  ## c(0, 4) is normalised to multipliers {0, 2} with mean 1
  p <- power_twfe(df, "y", "x", "id", "t", reps = 60, het = c(0, 4))
  expect_match(p$extras$het$label, "\\{0, 2\\}")
})

test_that("het runs across the other designs (smoke)", {
  set.seed(304)
  iv <- sim_iv(n = 400)
  pv <- power_iv(iv, "y", "d", z = "z", reps = 150, het = c(0, 2))
  expect_true(is.finite(pv$mde["80%"]))
  expect_gte(pv$mde["80%"], pv$extras$het$mde_constant["80%"])
  rd <- power_rdd(sim_rdd(n = 800), "y", run = "x", cutoff = 0, bandwidth = 0.5,
                  reps = 150, bw_reps = 40, het = c(0, 2))
  expect_true(is.finite(rd$mde["80%"]))
  ev <- sim_event(n_units = 60, n_time = 8, adopt = 5)
  pe <- power_event(ev, "y", "id", "t", treat = "d", reps = 150, het = c(0, 2))
  expect_true(is.finite(pe$mde["80%"]))
  expect_gte(pe$mde["80%"], pe$extras$het$mde_constant["80%"])
  ## hypothetical-instrument mode
  iv$cand <- rnorm(nrow(iv))
  ph <- power_iv(iv, "y", "d", L = "cand", rho = c(0.3, 0.5), reps = 120,
                 het = c(0, 2))
  expect_true(all(is.finite(ph$mde[["80%"]])))
})
