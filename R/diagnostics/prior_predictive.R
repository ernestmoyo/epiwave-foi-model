# ==============================================================================
# Prior predictive checks for the Stage 2 GP model
# ==============================================================================
# Punam's suggestion (2026-03-30 meeting): before fitting, draw the parameters
# from their priors, push them through the model, and check that the priors imply
# plausible data — and, critically, whether the spatial lengthscale prior lives in
# a region the data can actually identify.
#
# This is pure forward simulation from the priors. No greta / no MCMC.
#
# Usage:
#   source("R/epiwave-foi-model.R")          # Stage 1 + helpers (no greta needed)
#   source("R/diagnostics/prior_predictive.R")
#   ppc <- run_prior_predictive(n_draws = 500)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

# --- Priors, exactly as declared in fit_epiwave_gp() ---------------------------
# alpha  ~ normal(0, 1)
# gamma  ~ normal(0.1, 0.05) truncated to (0.001, Inf)
# sigma2 ~ lognormal(-0.5, 0.5)
# phi    ~ lognormal(-0.7, 0.5)   (median ~0.5; domain-aware, identifiable range)
# theta  ~ uniform(0, 1)
draw_priors <- function(n) {
  rtnorm <- function(n, m, s, lo) { x <- rnorm(n, m, s); x[x < lo] <- lo; x }
  data.frame(
    alpha  = rnorm(n, 0, 1),
    gamma  = rtnorm(n, 0.1, 0.05, 0.001),
    sigma2 = rlnorm(n, -0.5, 0.5),
    phi    = rlnorm(n, -0.7, 0.5),
    theta  = runif(n, 0, 1)
  )
}

# --- Matern 5/2 correlation as a function of distance and lengthscale ----------
matern52_corr <- function(d, phi) {
  r <- sqrt(5) * d / phi
  (1 + r + r^2 / 3) * exp(-r)
}

#' Run the prior predictive check.
#'
#' @param n_draws Number of prior draws
#' @param n_sites,n_times Grid size (matches the study)
#' @param seed Reproducibility seed
#' @param outdir Where to write the plots
#' @return A list with the per-draw summary table and the fixed mechanistic pieces
#' @export
run_prior_predictive <- function(n_draws = 500, n_sites = 10, n_times = 48,
                                 seed = 2026,
                                 outdir = "outputs/prior_predictive") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)

  # Fixed mechanistic scaffold (I*, coords, convolution, population) from one call
  # to the data generator. We only reuse its deterministic pieces; the GP field is
  # drawn fresh from the priors below.
  d <- simulate_epiwave_data(n_sites = n_sites, n_times = n_times, seed = 1)
  I_star_t   <- t(pmax(d$I_star, 1e-6))      # [n_sites x n_times]
  N_t        <- t(d$pop_matrix)              # [n_sites x n_times]
  coords     <- d$spatial_coords_norm        # [n_sites x 2], on [0,1]
  conv       <- d$conv_matrix                # [n_times x n_times]
  nt         <- ncol(I_star_t)
  max_dist   <- max(dist(coords))            # largest inter-site distance on [0,1]
  true_phi   <- d$true_params$gp_phi         # the lengthscale the study tries to recover

  pri <- draw_priors(n_draws)

  # min pairwise correlation implied by phi = correlation at the largest distance.
  # If this is near 1, every site is ~identically correlated and phi is
  # unidentifiable no matter how long we sample.
  pri$min_corr <- matern52_corr(max_dist, pri$phi)

  # Push each prior draw through the model and record what data it implies
  summ <- lapply(seq_len(n_draws), function(i) {
    p <- pri[i, ]
    eps <- simulate_gp_residuals(coords, nt, sigma = sqrt(p$sigma2),
                                 phi = p$phi, rho = p$theta)      # [n_sites x n_times]
    I   <- exp(p$alpha + log(I_star_t) + eps)                    # latent incidence
    prev <- 1 - exp(-(I %*% conv))                               # observed-scale prevalence
    exp_cases <- p$gamma * I * N_t                               # expected case counts
    data.frame(
      draw = i, phi = p$phi, sigma2 = p$sigma2, min_corr = p$min_corr,
      prev_median = median(prev), prev_max = max(prev),
      prev_frac_saturated = mean(prev > 0.99),
      cases_median = median(exp_cases), cases_max = max(exp_cases)
    )
  })
  summ <- do.call(rbind, summ)

  # ---- Diagnostic numbers ----
  frac_flat <- mean(pri$min_corr > 0.8)
  frac_sat  <- mean(summ$prev_frac_saturated > 0.10)
  cat(sprintf("\nPrior draws where ALL sites are >80%% correlated (phi unidentifiable): %.0f%%\n",
              100 * frac_flat))
  cat(sprintf("Prior draws where >10%% of prevalence cells saturate at ~1: %.0f%%\n",
              100 * frac_sat))
  cat(sprintf("Median prior-predictive prevalence: %.3f  |  median expected cases: %.1f\n",
              median(summ$prev_median), median(summ$cases_median)))

  # ---- Plots ----
  # 1. phi prior -> minimum spatial correlation (the identifiability diagnostic)
  p1 <- ggplot(pri, aes(min_corr)) +
    geom_histogram(bins = 40, fill = "#2563eb", colour = "white") +
    geom_vline(xintercept = 0.8, linetype = "dashed", colour = "#b91c1c") +
    annotate("text", x = 0.8, y = Inf, label = "  >0.8: sites indistinguishable",
             hjust = 0, vjust = 1.5, size = 3.4, colour = "#b91c1c") +
    labs(title = sprintf("%.0f%% of the lengthscale prior sits in the unidentifiable region",
                         100 * mean(pri$min_corr > 0.8)),
         subtitle = sprintf("Min correlation between the two farthest sites, over %d prior draws of phi", n_draws),
         x = "Minimum pairwise spatial correlation implied by phi", y = "Prior draws") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "ppc_phi_correlation.png"), p1, width = 7.6, height = 4.3, dpi = 150)

  # 2. Prior-predictive prevalence (is the data plausible?)
  p2 <- ggplot(summ, aes(prev_median)) +
    geom_histogram(bins = 40, fill = "#b91c1c", colour = "white") +
    annotate("rect", xmin = 0, xmax = 0.4, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#16a34a") +
    annotate("text", x = 0.2, y = Inf, label = "plausible range", vjust = 1.5, size = 3.4, colour = "#16a34a") +
    labs(title = sprintf("The priors imply implausibly high prevalence (median %.0f%%)",
                         100 * median(summ$prev_median)),
         subtitle = "Median prevalence per prior draw piles up near 1.0 — the incidence scale is too high",
         x = "Median prevalence implied by a prior draw", y = "Prior draws") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "ppc_prevalence.png"), p2, width = 7.6, height = 4.3, dpi = 150)

  # 3. What lengthscale WOULD be identifiable — correlation vs phi
  phi_grid <- seq(0.05, 4, length.out = 200)
  curve_df <- data.frame(phi = phi_grid, min_corr = matern52_corr(max_dist, phi_grid))
  phi_star <- approx(curve_df$min_corr, curve_df$phi, xout = 0.5)$y  # phi giving 0.5 corr
  p3 <- ggplot(curve_df, aes(phi, min_corr)) +
    geom_line(linewidth = 1, colour = "#2563eb") +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#b91c1c") +
    geom_vline(xintercept = true_phi, linetype = "dotted", colour = "#16a34a") +
    annotate("text", x = true_phi, y = 0.6, label = sprintf("simulation now uses phi = %.1f", true_phi),
             hjust = -0.08, size = 3.4, colour = "#16a34a") +
    labs(title = "Only short lengthscales are identifiable on this grid",
         subtitle = sprintf("Correlation at the largest site distance; identifiable structure needs phi below ~%.1f", phi_star),
         x = "Lengthscale phi", y = "Min pairwise correlation") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "ppc_identifiable_range.png"), p3, width = 7.6, height = 4.3, dpi = 150)

  cat(sprintf("phi giving 0.5 correlation at max distance: %.2f (prior median %.2f, sim truth %.2f)\n",
              phi_star, median(pri$phi), true_phi))
  cat("Saved plots to", outdir, "\n")

  invisible(list(priors = pri, summary = summ, max_dist = max_dist, phi_star = phi_star))
}
