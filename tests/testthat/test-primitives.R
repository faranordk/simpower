test_that("two-way within estimator matches lm with factor FE", {
  set.seed(1)
  df <- sim_panel(n_units = 40, n_time = 8, beta = 0.5)
  est <- simpower:::.fit_twfe(df$y, df$x,
                             as.integer(factor(df$id)),
                             as.integer(factor(df$t)))
  ref <- stats::lm(y ~ x + factor(id) + factor(t), data = df)
  expect_equal(unname(est["b"]), unname(coef(ref)["x"]), tolerance = 1e-6)
})

test_that("2SLS matches AER::ivreg coefficient", {
  skip_if_not_installed("AER")
  set.seed(2)
  df <- sim_iv(n = 500, strength = 0.6, beta = 0.4)
  est <- simpower:::.fit_tsls(df$y, df$d, matrix(df$z, ncol = 1))
  ref <- AER::ivreg(y ~ d | z, data = df)
  expect_equal(unname(est["b"]), unname(coef(ref)["d"]), tolerance = 1e-8)
})

test_that("power curve is monotone and starts near alpha", {
  set.seed(3)
  bs <- rnorm(2000, 0, 1); ses <- rep(1, 2000)
  grid <- seq(0, 5, by = 0.05)
  p <- power_from_null(bs, ses, grid, alpha = 0.05, alternative = "greater")
  expect_true(all(diff(p) >= -1e-8))          # non-decreasing
  expect_lt(abs(p[1] - 0.05), 0.02)            # power at 0 ~ alpha
  expect_gt(p[length(p)], 0.99)                # saturates
})

test_that("one-sided 'less' injects effects in the negative direction", {
  set.seed(4)
  bs <- rnorm(2000, 0, 1); ses <- rep(1, 2000)
  grid <- seq(0, 5, by = 0.05)
  p <- power_from_null(bs, ses, grid, alpha = 0.05, alternative = "less")
  expect_true(all(diff(p) >= -1e-8))           # power rises with effect magnitude
  expect_lt(abs(p[1] - 0.05), 0.02)
  expect_gt(p[length(p)], 0.99)
})

test_that("uniform-kernel RDD fit equals OLS within the window", {
  set.seed(5)
  df <- sim_rdd(n = 1200, jump = 0.4)
  inb <- abs(df$x) < 0.5
  est <- simpower:::.fit_rdd(df$y[inb], df$x[inb], as.numeric(df$x[inb] >= 0),
                            w = rep(1, sum(inb)))
  ref <- stats::lm(y ~ I(as.numeric(x >= 0)) + x + I(as.numeric(x >= 0) * x),
                   data = df[inb, ])
  expect_equal(unname(est["b"]), unname(coef(ref)[2]), tolerance = 1e-7)
})

test_that("mde_from_curve finds the crossing", {
  grid <- seq(0, 1, by = 0.1)
  power <- seq(0, 1, by = 0.1)
  expect_equal(mde_from_curve(grid, power, 0.8), 0.8)
})
