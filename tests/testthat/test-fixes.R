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
  sm <- simpower:::.subsample_units(uid, treated_ids = 1:5,
                                    n_units = NULL, n_treated = NULL)
  expect_true(all(sm))
  expect_length(sm, length(uid))
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
