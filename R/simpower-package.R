#' simpower: Simulation-Based Statistical Power for Causal Designs
#'
#' `simpower` computes statistical power and the minimum detectable effect
#' (MDE) for four causal-inference designs --- two-way fixed effects
#' ([power_twfe()]), event study ([power_event()]), instrumental variables
#' ([power_iv()]), and regression discontinuity ([power_rdd()]) --- directly
#' from a user's own data, following the Monte-Carlo recombination method of
#' Doucette (2025).
#'
#' The idea is the same across designs. The outcome is repeatedly *recombined*
#' so that any true relationship to the treatment is broken (the empirical
#' null). The design's own estimator and standard errors are computed on each
#' recombined sample and stored. Because every estimator here is linear in the
#' outcome, a growing true effect can be injected *analytically*: adding a true
#' effect \eqn{\tau} shifts the estimated coefficient by exactly \eqn{\tau} and
#' leaves the standard error unchanged, so the whole power curve is a cheap
#' post-processing step over a single set of simulations.
#'
#' This preserves the real features of the data that determine power ---
#' autocorrelation, clustering, unbalanced panels, instrument strength, and the
#' shape of the running variable --- rather than assuming them away with a
#' closed-form formula.
#'
#' @section Main functions:
#' \describe{
#'   \item{[power_twfe()]}{Two-way (unit + time) fixed effects.}
#'   \item{[power_event()]}{Event study with a choice of effect shape and a
#'     per-horizon MDE graph.}
#'   \item{[power_iv()]}{Instrumental variables / two-stage least squares.}
#'   \item{[power_rdd()]}{Sharp regression discontinuity, with an
#'     MDE-versus-bandwidth graph.}
#' }
#'
#' @references
#' Doucette, M. (2025). What Can We Learn About the Effects of Democracy Using
#' Cross-National Data? American Political Science Review, 119(3).
#'
#' Doucette, M. (2026). simpower: Simulation-Based Statistical Power for Panel
#' and Causal Designs. arXiv.
#'
#' @keywords internal
"_PACKAGE"

## quiet R CMD check for ggplot2 aes() variables used via .data is unnecessary,
## but we declare NULL bindings for any bare-symbol columns just in case.
utils::globalVariables(c("effect", "power", "horizon", "mde", "mde_peak",
                         "bandwidth", "lower", "upper", "kind", "xval", "rho"))

`%||%` <- function(a, b) if (is.null(a)) b else a
