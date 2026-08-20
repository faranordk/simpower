#' @export
print.simpower <- function(x, ...) {
  designlab <- c(twfe = "Two-way fixed effects", event = "Event study",
                 iv = "Instrumental variables (2SLS)",
                 rdd = "Regression discontinuity")[x$design]
  cat(sprintf("<simpower> %s\n", designlab))
  alt <- c(greater = "one-sided (greater)", less = "one-sided (less)",
           two.sided = "two-sided")[x$alternative]
  cat(sprintf("  test: %s, alpha = %.3g, reps = %d\n", alt, x$alpha, x$reps))
  sdy <- x$extras$sd_y
  if (is.data.frame(x$mde)) {
    ## IV hypothetical-instrument mode: one MDE row per assumed cor(z, x).
    ## The pooled "mean se" is meaningless here (dominated by the weakest
    ## instrument's occasional huge draws), so report a per-rho median SE.
    cat(sprintf("  hypothetical instrument: L = '%s', curves by cor(z, x)\n",
                x$extras$L))
    tab <- x$mde
    tab$med_se <- apply(x$null$ses, 2L,
                        function(s) stats::median(s[is.finite(s) & s > 0]))
    tab$implied_F <- unname(x$extras$first_stage_F)
    if (!is.null(sdy) && is.finite(sdy) && sdy > 0) {
      for (nm in setdiff(names(x$mde), "rho")) {
        tab[[paste0(nm, " (SD)")]] <- x$mde[[nm]] / sdy
      }
    }
    cat("  MDE by instrument strength:\n")
    print(format(tab, digits = 3), row.names = FALSE)
  } else {
    cat(sprintf("  null: mean beta = %.4f, mean se = %.4f\n",
                x$null$b_mean, x$null$se_mean))
    mline <- paste(sprintf("%s power = %.3g", names(x$mde), x$mde), collapse = "   ")
    cat(sprintf("  MDE:  %s\n", mline))
    if (!is.null(x$extras$het) && !is.null(x$extras$het$mde_constant)) {
      cat(sprintf("  heterogeneous effects: %s\n", x$extras$het$label))
      pc <- x$extras$het$power_ceiling
      if (!is.null(pc) && is.finite(pc) && pc < 0.995) {
        cat(sprintf("  power ceiling under this heterogeneity: %.2f\n", pc))
      }
      mc <- x$extras$het$mde_constant
      cat(sprintf("  MDE if constant (for comparison):  %s\n",
                  paste(sprintf("%s power = %.3g", names(mc), mc), collapse = "   ")))
    }
    if (!is.null(sdy) && is.finite(sdy) && sdy > 0) {
      sline <- paste(sprintf("%s power = %.3g", names(x$mde), x$mde / sdy), collapse = "   ")
      cat(sprintf("  MDE (SD units):  %s   [within-SD of Y = %.4g]\n", sline, sdy))
    }
  }
  if (x$design == "event") {
    cat(sprintf("  effect shape: %s   (MDE in units of the peak effect)\n",
                x$extras$shape))
  }
  if (x$design == "rdd") {
    cat(sprintf("  bandwidth: %.4g  (%s)\n", x$extras$h, x$extras$bw_source))
  }
  if (x$design == "iv" && all(is.finite(x$extras$estimate))) {
    cat(sprintf("  first-stage F: %.1f\n", x$extras$first_stage_F))
  }
  if (x$design == "twfe" && !is.null(x$extras$estimate) &&
      all(is.finite(x$extras$estimate))) {
    cat(sprintf("  real-data estimate: b = %.4f (se = %.4f, %s)\n",
                x$extras$estimate["b"], x$extras$estimate["se"], x$extras$estimator))
  }
  cat("  use summary() for the sample-composition check and plot() for figures.\n")
  invisible(x)
}

#' Summarise a simpower object
#' @param object A `simpower` object.
#' @param ... Unused.
#' @export
summary.simpower <- function(object, ...) {
  print(object)
  cat("\nSample check (target vs simulation average):\n")
  sc <- object$sample_check
  wid <- max(nchar(sc$metric)) + 2L        # pad to the longest label
  for (i in seq_len(nrow(sc))) {
    cat(sprintf("  %-*s target = %-8.1f sim avg = %.1f\n",
                wid, sc$metric[i], sc$target[i], sc$sim_avg[i]))
  }
  if (object$design == "event") {
    cat("\nPer-horizon MDE (at ",
        attr(object$extras$horizon, "target") * 100, "% power):\n", sep = "")
    print(format(object$extras$horizon, digits = 3), row.names = FALSE)
  }
  if (object$design == "rdd") {
    cat("\nMDE by bandwidth (at ",
        attr(object$extras$bandwidth, "target") * 100, "% power):\n", sep = "")
    print(format(object$extras$bandwidth, digits = 3), row.names = FALSE)
  }
  invisible(object)
}

#' Plot a simpower object
#'
#' @param x A `simpower` object.
#' @param type Which figure(s) to draw. For event studies: `"overall"`,
#'   `"horizon"`, or `"both"` (default). For RDD: `"overall"`, `"bandwidth"`,
#'   or `"both"` (default). Ignored for TWFE/IV (single curve).
#' @param ... Unused.
#' @return A ggplot (or patchwork) object.
#' @export
plot.simpower <- function(x, type = c("both", "overall", "horizon", "bandwidth"),
                          ...) {
  type <- match.arg(type)
  if (x$design %in% c("twfe", "iv")) return(.gg_curve(x))
  if (x$design == "event") {
    if (type == "overall")  return(.gg_curve(x))
    if (type == "horizon")  return(.gg_horizon(x))
    ## Combine via wrap_plots(), not `p1 + p2`: the `+` method for two ggplots
    ## is only registered once the patchwork namespace is loaded, and in
    ## `p1 + p2 + patchwork::plot_annotation(...)` the first `+` is evaluated
    ## *before* `patchwork::` triggers that load -- so plot() would error in a
    ## fresh session unless the user had attached patchwork themselves.
    return(patchwork::wrap_plots(.gg_curve(x), .gg_horizon(x)) +
             patchwork::plot_annotation(title = "Event-study power"))
  }
  if (x$design == "rdd") {
    if (type == "overall")   return(.gg_curve(x))
    if (type == "bandwidth") return(.gg_bandwidth(x))
    return(patchwork::wrap_plots(.gg_curve(x), .gg_bandwidth(x)) +
             patchwork::plot_annotation(title = "RDD power"))
  }
  .gg_curve(x)
}

.pw_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

.gg_curve <- function(x) {
  df <- x$power
  if (!is.null(x$extras$rho)) {
    ## IV hypothetical-instrument mode: one curve per assumed cor(z, x)
    th <- as.numeric(sub("%", "", setdiff(names(x$mde), "rho")[1L])) / 100
    df$rho <- factor(df$rho)
    return(
      ggplot2::ggplot(df, ggplot2::aes(effect, power, colour = rho)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_hline(yintercept = th, linetype = 2, colour = "grey50") +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(title = "IV power by instrument strength",
                      subtitle = sprintf("hypothetical instrument built from '%s'",
                                         x$extras$L),
                      x = "True effect size", y = "Power (share significant)",
                      colour = "cor(z, x)") +
        .pw_theme()
    )
  }
  m80 <- x$mde[1]; th80 <- as.numeric(sub("%", "", names(x$mde)[1])) / 100
  ttl <- switch(x$design,
                twfe  = "TWFE power curve",
                iv    = "IV (2SLS) power curve",
                event = sprintf("Overall power (%s effect)", x$extras$shape),
                rdd   = sprintf("RDD power (bandwidth = %.3g)", x$extras$h))
  g <- ggplot2::ggplot(df, ggplot2::aes(effect, power)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = th80, linetype = 2, colour = "grey50") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = ttl,
                  x = "True effect size", y = "Power (share significant)") +
    .pw_theme()
  if (is.finite(m80)) {
    g <- g + ggplot2::geom_vline(xintercept = m80, linetype = 3, colour = "steelblue") +
      ggplot2::annotate("point", x = m80, y = th80, colour = "steelblue", size = 2)
  }
  g
}

.gg_horizon <- function(x) {
  hz    <- x$extras$horizon
  th    <- attr(hz, "target")
  shape <- x$extras$shape
  ## Take the assumed effect shape into account: plot the per-horizon MDE in
  ## units of the *peak* effect (mde / weight), matching the overall MDE. For a
  ## "constant" shape the weights are all 1, so this equals the raw per-horizon
  ## MDE and the plot is unchanged.
  ylab <- if (identical(shape, "constant"))
            "Minimum detectable effect"
          else "Per-horizon MDE (peak-effect units)"
  ggplot2::ggplot(hz, ggplot2::aes(factor(horizon), mde_peak)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.85) +
    ggplot2::labs(title = sprintf("Per-horizon MDE (%.0f%% power, %s shape)",
                                  th * 100, shape),
                  x = "Event-time horizon", y = ylab) +
    .pw_theme()
}

.gg_bandwidth <- function(x) {
  bw <- x$extras$bandwidth
  th <- attr(bw, "target")
  ggplot2::ggplot(bw, ggplot2::aes(bandwidth, mde)) +
    ggplot2::geom_line(linewidth = 1, colour = "steelblue") +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_vline(xintercept = x$extras$h, linetype = 3) +
    ggplot2::labs(title = sprintf("MDE vs bandwidth (%.0f%% power)", th * 100),
                  x = "Bandwidth", y = "Minimum detectable effect") +
    .pw_theme()
}

#' Extract the minimum detectable effect(s)
#'
#' @param object A `simpower` object.
#' @param target Optional power level(s) at which to read off the MDE. When
#'   `NULL` (default) the targets computed at fit time are returned.
#' @param standardized If `TRUE`, divide by the within-unit SD of the outcome to
#'   express the MDE in standard-deviation units.
#' @param ... Unused.
#' @return Named numeric vector of MDEs. For IV objects fitted in
#'   hypothetical-instrument mode (`power_iv(L = ...)`), a data frame with one
#'   row per `rho` instead.
#' @export
mde <- function(object, target = NULL, standardized = FALSE, ...) UseMethod("mde")

#' @export
mde.simpower <- function(object, target = NULL, standardized = FALSE, ...) {
  sdy <- object$extras$sd_y
  sdy_ok <- !is.null(sdy) && is.finite(sdy) && sdy > 0
  if (standardized && !sdy_ok)
    stop("No usable within-SD of the outcome is stored.", call. = FALSE)
  if (is.data.frame(object$mde)) {
    ## IV hypothetical-instrument mode: return a data frame, one row per rho
    m <- if (is.null(target)) {
      object$mde
    } else {
      rhos <- object$extras$rho
      cbind(data.frame(rho = rhos), do.call(rbind, lapply(rhos, function(r) {
        p <- object$power$power[object$power$rho == r]
        as.data.frame(as.list(stats::setNames(
          vapply(target, function(th) mde_from_curve(object$grid, p, th), numeric(1)),
          paste0(target * 100, "%"))), check.names = FALSE)
      })))
    }
    if (standardized) {
      for (nm in setdiff(names(m), "rho")) m[[nm]] <- m[[nm]] / sdy
    }
    return(m)
  }
  m <- if (is.null(target)) {
    object$mde
  } else {
    stats::setNames(vapply(target,
      function(th) mde_from_curve(object$grid, object$power$power, th),
      numeric(1)), paste0(target * 100, "%"))
  }
  if (standardized) m <- m / sdy
  m
}

#' Power at a given true effect size
#'
#' Reads the simulated power curve at arbitrary effect value(s) by linear
#' interpolation.
#'
#' @param object A `simpower` object.
#' @param effect Numeric vector of true effect sizes (magnitudes in the
#'   hypothesised direction).
#' @param ... Unused.
#' @return Numeric vector of power values. For IV objects fitted in
#'   hypothetical-instrument mode, a data frame with an `effect` column and one
#'   power column per `rho`.
#' @examples
#' set.seed(1)
#' pw <- power_twfe(sim_panel(50, 8), "y", "x", "id", "t", reps = 120)
#' power_at(pw, c(0.1, 0.2, 0.3))
#' @export
power_at <- function(object, effect, ...) UseMethod("power_at")

#' @export
power_at.simpower <- function(object, effect, ...) {
  if (!is.null(object$extras$rho)) {
    ## IV hypothetical-instrument mode: one power column per rho
    out <- data.frame(effect = effect)
    for (r in object$extras$rho) {
      sub <- object$power[object$power$rho == r, ]
      out[[paste0("r", format(r, trim = TRUE, drop0trailing = TRUE))]] <-
        stats::approx(sub$effect, sub$power, xout = effect, rule = 2)$y
    }
    return(out)
  }
  stats::approx(object$power$effect, object$power$power, xout = effect, rule = 2)$y
}

#' Extract simpower tables as a data frame
#'
#' @param x A `simpower` object.
#' @param row.names,optional Ignored; present for S3 consistency with
#'   [base::as.data.frame()].
#' @param what Which table: `"power"` (the power curve, default), `"horizon"`
#'   (event study), or `"bandwidth"` (RDD).
#' @param ... Unused.
#' @return A data frame.
#' @export
as.data.frame.simpower <- function(x, row.names = NULL, optional = FALSE,
                                   what = c("power", "horizon", "bandwidth"), ...) {
  what <- match.arg(what)
  if (what == "power") return(x$power)
  if (what == "horizon") {
    if (x$design != "event") stop("`horizon` is only available for event studies.", call. = FALSE)
    return(as.data.frame(unclass(x$extras$horizon)))
  }
  if (x$design != "rdd") stop("`bandwidth` is only available for RDD.", call. = FALSE)
  as.data.frame(unclass(x$extras$bandwidth))
}
