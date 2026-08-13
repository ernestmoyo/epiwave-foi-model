# ==============================================================================
# Sensitivity of I* to the fixed entomological parameters
# ==============================================================================
# David's suggestion (Melbourne, 12 Aug 2026): the entomological inputs m, a, g
# are treated as fixed, but they are posterior summaries from the Vector Atlas
# models, not known constants. Rather than propagating full uncertainty (which
# would reintroduce the computational cost the two-stage design avoids), run a
# cheap sensitivity: solve Stage 1 with a lower and an upper value in place of
# the central one, and report how much I* moves.
#
# The scenarios below scale the central value by +/- 20 percent as a stand-in
# for the 25th and 75th posterior percentiles. When the Vector Atlas surfaces
# arrive, replace the scaling with the actual posterior quantiles.
#
# Stage 1 only. No MCMC, runs in seconds.
#
# Usage:
#   source("R/epiwave-foi-model.R")
#   source("R/diagnostics/sensitivity_fixed_params.R")
#   sens <- run_fixed_param_sensitivity()
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

#' Solve Stage 1 once for a given scaling of the entomological parameters.
#' @keywords internal
.solve_istar <- function(scale_m = 1, scale_a = 1, scale_g = 1,
                         n_sites = 10, n_times = 48, seed = 123) {
  times     <- seq(0, n_times * 30, by = 30)
  locations <- paste0("Site_", sprintf("%02d", 1:n_sites))
  m <- get_fixed_m(times, locations, baseline_m = 2.0 * scale_m,
                   seasonal_amplitude = 0.6)
  a <- get_fixed_a(times, locations, baseline_a = 0.3 * scale_a)
  g <- get_fixed_g(times, locations, baseline_g = (1/10) * scale_g)
  ode <- solve_ross_macdonald_multi_site(m, a, g, times = times,
                                         b = 0.8, c = 0.8, r = 1/7)
  compute_mechanistic_prediction(m, a, 0.8, ode$z)
}

#' Run the one-at-a-time sensitivity of I* to m, a and g.
#'
#' @param spread Proportional spread around the central value (default 0.20,
#'   standing in for the 25th/75th posterior percentiles)
#' @param n_sites,n_times Grid size (matches the study)
#' @param site_idx Site shown in the time-series plot (default 1)
#' @param outdir Where to write the plot and table
#' @return List with the summary table and the per-scenario I* matrices
#' @export
run_fixed_param_sensitivity <- function(spread = 0.20,
                                        n_sites = 10, n_times = 48,
                                        site_idx = 1,
                                        outdir = "outputs/sensitivity") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  lo <- 1 - spread; hi <- 1 + spread

  scenarios <- list(
    central  = c(1, 1, 1),
    m_low    = c(lo, 1, 1),  m_high  = c(hi, 1, 1),
    a_low    = c(1, lo, 1),  a_high  = c(1, hi, 1),
    g_low    = c(1, 1, lo),  g_high  = c(1, 1, hi)
  )

  istars <- lapply(scenarios, function(s)
    .solve_istar(s[1], s[2], s[3], n_sites = n_sites, n_times = n_times))

  central_mean <- mean(istars$central)
  summary_tab <- do.call(rbind, lapply(names(scenarios), function(nm) {
    data.frame(scenario = nm,
               mean_istar = mean(istars[[nm]]),
               pct_change = 100 * (mean(istars[[nm]]) - central_mean) / central_mean)
  }))
  cat("\nSensitivity of mean I* to the fixed entomological parameters",
      sprintf("(+/- %.0f%% around the central value):\n", 100 * spread))
  print(summary_tab, row.names = FALSE, digits = 3)

  # Time-series plot at one site: central line with low/high band per parameter
  mk_df <- function(par) {
    data.frame(
      month   = seq_len(nrow(istars$central)),
      central = istars$central[, site_idx],
      low     = istars[[paste0(par, "_low")]][, site_idx],
      high    = istars[[paste0(par, "_high")]][, site_idx],
      param   = par)
  }
  band_df <- do.call(rbind, lapply(c("m", "a", "g"), mk_df))
  band_df$param <- factor(band_df$param, levels = c("m", "a", "g"),
                          labels = c("m (abundance)", "a (biting rate)", "g (mortality)"))

  p <- ggplot(band_df, aes(x = month)) +
    geom_ribbon(aes(ymin = pmin(low, high), ymax = pmax(low, high)),
                fill = "#2E75B6", alpha = 0.25) +
    geom_line(aes(y = central), colour = "#1F3A5F", linewidth = 0.9) +
    facet_wrap(~ param, ncol = 1) +
    labs(title = sprintf("How much I* moves when a fixed parameter shifts by %.0f%%", 100 * spread),
         subtitle = sprintf("Site %d. The band spans the low and high scenarios; the line is the central value.", site_idx),
         x = "Month", y = "I* (mechanistic incidence rate)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
  ggsave(file.path(outdir, "sensitivity_istar.png"), p, width = 7.6, height = 6.4, dpi = 150)

  write.csv(summary_tab, file.path(outdir, "sensitivity_summary.csv"), row.names = FALSE)
  cat("Saved plot and table to", outdir, "\n")
  invisible(list(summary = summary_tab, istars = istars))
}
