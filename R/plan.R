#' Trace the MDE across sample sizes (design planning)
#'
#' Re-runs a fitted [simpower] analysis over a grid of `n_units` and/or
#' `n_treated` values and returns how the minimum detectable effect changes with
#' sample size --- the core question in ex-ante design planning ("how many units
#' do I need to detect an effect of size X?").
#'
#' @param object A fitted `simpower` object (from any of the `power_*`
#'   functions). Its stored call is re-evaluated with the sample-size arguments
#'   replaced.
#' @param n_units,n_treated Vectors of sample sizes to sweep. Supply either or
#'   both; if both are given they are recycled to a common length and varied in
#'   lockstep (pairs of total/treated counts).
#' @param target Power level at which to read the MDE.
#' @param standardized If `TRUE`, report the MDE in within-SD-of-outcome units.
#' @param reps Optional override for the number of simulations per point (useful
#'   to keep the sweep fast).
#' @param plot If `TRUE`, the returned object prints as an MDE-vs-sample-size
#'   plot.
#' @param envir Environment in which to re-evaluate the original call (defaults
#'   to the caller); must be able to see the original `data`.
#'
#' @return An object of class `"simpower_plan"`: a list with `table` (a data
#'   frame of sample size and MDE) and `plot` (a ggplot).
#'
#' @examples
#' set.seed(1)
#' df <- sim_panel(120, 8)
#' pw <- power_twfe(df, "y", "x", "id", "t", reps = 100)
#' pl <- plan_mde(pw, n_units = c(30, 60, 90, 120), reps = 100)
#' pl$table
#' @export
plan_mde <- function(object, n_units = NULL, n_treated = NULL, target = 0.8,
                     standardized = FALSE, reps = NULL, plot = TRUE,
                     envir = parent.frame()) {
  stopifnot(inherits(object, "simpower"))
  if (!is.null(object$extras$rho)) {
    stop(paste0("`plan_mde()` does not support hypothetical-instrument IV objects ",
                "(one MDE per `rho`); call `power_iv()` with different `n_units` ",
                "values directly."), call. = FALSE)
  }
  if (is.null(n_units) && is.null(n_treated)) {
    stop("Supply `n_units` and/or `n_treated` as a vector of sample sizes.",
         call. = FALSE)
  }
  cl <- object$call
  L  <- max(length(n_units), length(n_treated))
  nu <- if (is.null(n_units))   vector("list", L) else as.list(rep_len(n_units, L))
  nt <- if (is.null(n_treated)) vector("list", L) else as.list(rep_len(n_treated, L))

  ## Prefer the data stored on the object at fit time; only fall back to
  ## re-evaluating the original `data` symbol in `envir` for objects created by
  ## older versions that did not store it.
  have_data <- !is.null(object$data)

  rows <- lapply(seq_len(L), function(i) {
    cl2 <- cl
    if (have_data)          cl2$data      <- object$data
    if (!is.null(nu[[i]])) cl2$n_units   <- nu[[i]]
    if (!is.null(nt[[i]])) cl2$n_treated <- nt[[i]]
    if (!is.null(reps))    cl2$reps      <- reps
    if (object$design == "rdd") {
      cl2$bw_grid <- object$extras$h      # skip the bandwidth sweep (unused here)
      cl2$bw_reps <- 2L                   # ...and don't pay for its reps
    }
    fit <- eval(cl2, envir)
    m  <- mde(fit, target = target, standardized = standardized)
    data.frame(n_units   = if (is.null(nu[[i]])) NA_real_ else nu[[i]],
               n_treated = if (is.null(nt[[i]])) NA_real_ else nt[[i]],
               mde       = unname(m[1]))
  })
  tab <- do.call(rbind, rows)

  xvar <- if (length(unique(stats::na.omit(tab$n_treated))) > 1 &&
              length(unique(stats::na.omit(tab$n_units))) <= 1) "n_treated" else "n_units"

  ## Build the ggplot only when it will be used -- so `plot = FALSE` does not
  ## require ggplot2 to be installed.
  gg <- NULL
  if (isTRUE(plot)) {
    xlab <- if (xvar == "n_treated") "Number of treated units" else "Number of units"
    ylab <- if (standardized) "MDE (SD units)" else "Minimum detectable effect"
    tab$xval <- tab[[xvar]]
    gg <- ggplot2::ggplot(tab, ggplot2::aes(xval, mde)) +
      ggplot2::geom_line(linewidth = 1, colour = "steelblue") +
      ggplot2::geom_point(colour = "steelblue") +
      ggplot2::labs(title = sprintf("MDE vs sample size (%.0f%% power, %s)",
                                    target * 100, object$design),
                    x = xlab, y = ylab) +
      ggplot2::scale_x_continuous(breaks = sort(unique(tab$xval))) +
      .pw_theme()
    tab$xval <- NULL
  }

  structure(list(table = tab, plot = gg, design = object$design,
                 target = target, standardized = standardized,
                 xvar = xvar, draw = isTRUE(plot)),
            class = "simpower_plan")
}

#' @export
print.simpower_plan <- function(x, ...) {
  cat(sprintf("<simpower_plan> MDE vs %s  (%s, %.0f%% power)\n",
              x$xvar, x$design, x$target * 100))
  print(format(x$table, digits = 3), row.names = FALSE)
  if (isTRUE(x$draw)) print(x$plot)
  invisible(x)
}

#' @export
plot.simpower_plan <- function(x, ...) x$plot
