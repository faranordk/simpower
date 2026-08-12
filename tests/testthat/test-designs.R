test_that("power_twfe returns a valid object with null centred at 0", {
  set.seed(10)
  df <- sim_panel(n_units = 50, n_time = 8)
  pw <- power_twfe(df, "y", "x", "id", "t", reps = 150)
  expect_s3_class(pw, "simpower")
  expect_equal(pw$design, "twfe")
  expect_lt(abs(pw$null$b_mean), 0.05)          # recombination imposes the null
  expect_true(is.finite(pw$mde["80%"]))
  expect_true(all(pw$power$power >= -1e-9 & pw$power$power <= 1 + 1e-9))
})

test_that("power_twfe responds to controls and subsampling", {
  set.seed(11)
  df <- sim_panel(n_units = 60, n_time = 8)
  df$z1 <- rnorm(nrow(df))
  pw <- power_twfe(df, "y", "x", "id", "t", controls = "z1",
                   n_units = 30, n_treated = 30, reps = 120)
  expect_lt(pw$sample_check$sim_avg[1], 60)     # subsampled to ~30 units
})

test_that("event study produces horizon table and shape matters", {
  set.seed(12)
  df <- sim_event(n_units = 80, n_time = 10, adopt = 6)
  pc <- power_event(df, "y", "id", "t", treat = "d", shape = "constant", reps = 120)
  pi <- power_event(df, "y", "id", "t", treat = "d", shape = "increasing", reps = 120)
  expect_equal(pc$design, "event")
  expect_true(nrow(pc$extras$horizon) >= 1)
  expect_true(all(c("horizon", "mde") %in% names(pc$extras$horizon)))
  # increasing shape puts less weight on early horizons -> different overall MDE
  expect_false(isTRUE(all.equal(unname(pc$mde["80%"]), unname(pi$mde["80%"]))))
})

test_that("power_iv computes first-stage F and centres the null", {
  set.seed(13)
  df <- sim_iv(n = 600, strength = 0.6)
  pw <- power_iv(df, "y", "d", "z", reps = 150)
  expect_equal(pw$design, "iv")
  expect_gt(pw$extras$first_stage_F, 10)        # strong instrument by construction
  expect_lt(abs(pw$null$b_mean), 0.1)
})

test_that("weak instruments give larger MDE than strong ones", {
  set.seed(14)
  strong <- power_iv(sim_iv(n = 600, strength = 0.8), "y", "d", "z", reps = 150)
  weak   <- power_iv(sim_iv(n = 600, strength = 0.15), "y", "d", "z", reps = 150)
  expect_gt(weak$mde["80%"], strong$mde["80%"])
})

test_that("standardized MDE, power_at, mde(target), and as.data.frame work", {
  set.seed(20)
  df <- sim_panel(60, 8)
  pw <- power_twfe(df, "y", "x", "id", "t", reps = 150)
  # standardized MDE = raw / within-SD
  expect_equal(unname(mde(pw, standardized = TRUE)["80%"]),
               unname(pw$mde["80%"]) / pw$extras$sd_y, tolerance = 1e-8)
  # arbitrary target
  m70 <- mde(pw, target = 0.7)
  expect_true(is.finite(m70) && m70 <= pw$mde["80%"])
  # power_at is increasing and bounded
  pa <- power_at(pw, c(0, pw$mde["80%"]))
  expect_lt(pa[1], pa[2]); expect_gte(pa[2], 0.79)
  # tidy accessor
  d <- as.data.frame(pw)
  expect_true(all(c("effect", "power") %in% names(d)))
})

test_that("a too-small emax warns instead of silently returning NA", {
  set.seed(30)
  df <- sim_panel(40, 6)
  expect_warning(
    pw <- power_twfe(df, "y", "x", "id", "t", emax = 0.001, step = 0.0005, reps = 100),
    "not reached within emax"
  )
  expect_true(is.na(pw$mde["80%"]))
})

test_that("plan_mde traces MDE down as units increase", {
  set.seed(21)
  df <- sim_panel(120, 8)
  pw <- power_twfe(df, "y", "x", "id", "t", reps = 100)
  pl <- plan_mde(pw, n_units = c(30, 60, 120), reps = 100, plot = FALSE)
  expect_s3_class(pl, "simpower_plan")
  expect_equal(nrow(pl$table), 3)
  expect_gt(pl$table$mde[1], pl$table$mde[3])   # more units -> smaller MDE
})

test_that("RDD returns a bandwidth sweep and MDE grows as bandwidth shrinks", {
  set.seed(15)
  df <- sim_rdd(n = 1500, jump = 0)
  pw <- suppressMessages(power_rdd(df, "y", "x", cutoff = 0, bandwidth = 0.5,
                                   reps = 150, bw_reps = 100))
  expect_equal(pw$design, "rdd")
  bt <- pw$extras$bandwidth
  expect_true(nrow(bt) >= 3)
  # smaller bandwidth -> fewer obs -> generally larger MDE (monotone-ish)
  expect_gt(bt$mde[which.min(bt$bandwidth)], bt$mde[which.max(bt$bandwidth)] * 0.9)
})
