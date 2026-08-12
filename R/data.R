## ---------------------------------------------------------------------------
## Small data simulators used in examples and tests. These generate toy data
## with realistic structure (unit effects, autocorrelation, endogeneity, a
## running-variable trend) so users can see each function work end-to-end.
## ---------------------------------------------------------------------------

#' Simulate a toy panel for TWFE examples
#'
#' @param n_units,n_time Number of units and time periods.
#' @param beta True coefficient on `x` (the recombination null ignores this;
#'   it only affects the observed data).
#' @param rho AR(1) parameter of the within-unit errors.
#' @param sd_fe,sd_e SD of the unit fixed effects and the idiosyncratic noise.
#' @return A data frame with columns `id`, `t`, `x`, `y`.
#' @export
sim_panel <- function(n_units = 100, n_time = 10, beta = 0.3, rho = 0.5,
                      sd_fe = 1, sd_e = 1) {
  id <- rep(seq_len(n_units), each = n_time)
  t  <- rep(seq_len(n_time), times = n_units)
  fe <- rep(stats::rnorm(n_units, 0, sd_fe), each = n_time)
  te <- rep(stats::rnorm(n_time, 0, 0.3), times = n_units)
  x  <- stats::rnorm(n_units * n_time)
  e  <- numeric(n_units * n_time)
  cf <- sqrt(1 - rho^2)
  for (i in seq_len(n_units)) {
    idx   <- ((i - 1L) * n_time + 1L):(i * n_time)   # rows for unit i (contiguous)
    innov <- stats::rnorm(n_time, 0, sd_e)
    xi <- innov; if (n_time >= 2L) xi[-1L] <- cf * innov[-1L]
    e[idx] <- as.numeric(stats::filter(xi, filter = rho, method = "recursive"))
  }
  y <- beta * x + fe + te + e
  data.frame(id = id, t = t, x = x, y = y)
}

#' Simulate a toy staggered-adoption panel for event-study examples
#'
#' @param n_units,n_time Number of units and time periods.
#' @param adopt Adoption period for treated units (half of units are treated).
#' @param rho AR(1) parameter of within-unit errors.
#' @return A data frame with columns `id`, `t`, `d` (absorbing 0/1 treatment), `y`.
#' @export
sim_event <- function(n_units = 100, n_time = 10, adopt = 6, rho = 0.5) {
  base <- sim_panel(n_units, n_time, beta = 0, rho = rho)
  treated_units <- sample(unique(base$id), floor(n_units / 2))
  base$d <- as.integer(base$id %in% treated_units & base$t >= adopt)
  base$x <- NULL
  base
}

#' Simulate toy cross-sectional data for IV examples
#'
#' @param n Number of observations.
#' @param strength First-stage strength (coefficient of `z` on `d`).
#' @param beta True structural effect of `d` on `y`.
#' @return A data frame with columns `y`, `d`, `z`.
#' @export
sim_iv <- function(n = 1000, strength = 0.5, beta = 0.4) {
  z <- stats::rnorm(n)
  conf <- stats::rnorm(n)                 # unobserved confounder
  d <- strength * z + conf + stats::rnorm(n)
  y <- beta * d + conf + stats::rnorm(n)  # confounder biases OLS, not IV
  data.frame(y = y, d = d, z = z)
}

#' Simulate toy data for RDD examples
#'
#' @param n Number of observations.
#' @param jump True discontinuity at the cutoff.
#' @param sd_e Noise SD.
#' @return A data frame with columns `x` (running variable, cutoff at 0), `y`.
#' @export
sim_rdd <- function(n = 2000, jump = 0.3, sd_e = 1) {
  x <- stats::runif(n, -1, 1)
  treat <- as.numeric(x >= 0)
  fx <- 0.8 * x + 0.5 * x^2            # smooth running-variable trend
  y  <- fx + jump * treat + stats::rnorm(n, 0, sd_e)
  data.frame(x = x, y = y)
}
