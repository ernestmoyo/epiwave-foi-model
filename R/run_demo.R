# ==============================================================================
# Worked example: does adding I* help?
# ==============================================================================
# Simulate one dataset from known truth, fit the Stage 2 model twice (with the
# I* offset, and with I* = 0, the standard geostatistical model), then check the
# fits against the truth. This is one replicate of the simulation study; the
# 50-replicate version is R/sim_estimation_harness.R.
#
# Run it line by line in RStudio from the project root, or with
# source("R/run_demo.R", echo = TRUE) so the plots print. The two fits take
# about 10-20 minutes on a laptop.
# ==============================================================================

Sys.setenv(RETICULATE_PYTHON = "C:/Users/ernes/AppData/Local/r-miniconda/envs/r-greta/python.exe")
source("R/greta_setup.R")
source("R/epiwave-foi-model.R")
library(ggplot2)

n_sites   <- 10
n_times   <- 48
n_samples <- 1000
warmup    <- 1000
chains    <- 2
site      <- 1     # the site shown in the time-series plots

# ---- 1. Simulate data from known truth (TRUE_PARAMS in the model file) -------
data <- simulate_epiwave_data(n_sites = n_sites, n_times = n_times)

# ---- 2. Fit the same model twice; only the offset differs --------------------
fit_model <- function(data, use_mechanistic) {
  model <- fit_epiwave_gp(
    observed_cases = data$observed_cases, I_star = data$I_star,
    N_pop = data$pop_matrix, spatial_coords = data$spatial_coords_norm,
    prev_data = data$prev_data, prev_conv_matrix = data$conv_matrix,
    use_mechanistic = use_mechanistic)
  draws <- mcmc(model, n_samples = n_samples, warmup = warmup, chains = chains)
  # posterior mean of the latent incidence surface, [n_sites x n_times]
  I_draws <- as.matrix(calculate(attr(model, "I_latent"), values = draws))
  list(draws = draws, I_hat = matrix(colMeans(I_draws), nrow = n_sites))
}

fit_with    <- fit_model(data, use_mechanistic = TRUE)
fit_without <- fit_model(data, use_mechanistic = FALSE)

# ---- 3. Score against the truth ------------------------------------------------
print(extract_posterior_summary(fit_with$draws))
print(coda::gelman.diag(fit_with$draws, multivariate = FALSE))  # R-hat: want < 1.05

scores <- rbind(
  with_I_star    = unlist(compute_performance_metrics(fit_with$I_hat, data$I_true_mat)),
  I_star_is_zero = unlist(compute_performance_metrics(fit_without$I_hat, data$I_true_mat)))
print(scores[, c("rmse", "mae", "relative_error")])   # lower is better


# ==============================================================================
# Plots
# ==============================================================================

theme_epiwave <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
months <- seq_along(data$times)

# The truth: how far the true incidence sits from I*. The gap is exp(epsilon).
ggplot(data.frame(month = months,
                  I_star = data$I_star[site, ],
                  I_true = data$I_true_mat[site, ]), aes(month)) +
  geom_line(aes(y = I_star, colour = "I* (mechanistic)"), linetype = "dashed") +
  geom_line(aes(y = I_true, colour = "true incidence")) +
  labs(title = sprintf("Site %d: I* and the true incidence", site),
       x = "Month", y = "Infections per person per day", colour = NULL) +
  theme_epiwave

# The true residual field epsilon, every site and month
ggplot(expand.grid(site = seq_len(n_sites), month = months) |>
         transform(epsilon = as.vector(data$epsilon_true_mat)),
       aes(month, factor(site), fill = epsilon)) +
  geom_tile() +
  scale_fill_gradient2(low = "#C00000", mid = "white", high = "#2E75B6") +
  labs(title = "True residual field epsilon", x = "Month", y = "Site") +
  theme_epiwave

# Alpha and gamma: a round cloud means separately identified; a diagonal
# band means only their product is identified.
posterior <- as.data.frame(as.matrix(fit_with$draws))
ggplot(posterior, aes(alpha, gamma_rr)) +
  geom_point(alpha = 0.15, colour = "#2E75B6") +
  geom_vline(xintercept = TRUE_PARAMS$alpha, linetype = "dashed") +
  geom_hline(yintercept = TRUE_PARAMS$reporting_rate, linetype = "dashed") +
  labs(title = "alpha and gamma (with I*)", x = "alpha", y = "gamma") +
  theme_epiwave

# Trace plots: chains that overlap and wander freely have mixed
trace <- do.call(rbind, lapply(seq_along(fit_with$draws), function(k) {
  d <- as.data.frame(as.matrix(fit_with$draws[[k]]))
  d$iteration <- seq_len(nrow(d)); d$chain <- factor(k); d
})) |>
  tidyr::pivot_longer(c(alpha, gamma_rr, sigma2, phi, theta),
                      names_to = "parameter", values_to = "value")
ggplot(trace, aes(iteration, value, colour = chain)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  labs(title = "Trace plots (with I*)", x = "Iteration after warmup", y = NULL) +
  theme_epiwave

# The headline: which fit recovers the true incidence?
ggplot(data.frame(month = months,
                  truth = data$I_true_mat[site, ],
                  with_I_star = fit_with$I_hat[site, ],
                  I_star_zero = fit_without$I_hat[site, ]), aes(month)) +
  geom_line(aes(y = truth, colour = "truth"), linewidth = 1) +
  geom_line(aes(y = with_I_star, colour = "with I*")) +
  geom_line(aes(y = I_star_zero, colour = "I* = 0"), linetype = "dashed") +
  scale_colour_manual(values = c("truth" = "black", "with I*" = "#2E75B6",
                                 "I* = 0" = "#C00000")) +
  labs(title = sprintf("Site %d: recovered incidence", site),
       x = "Month", y = "Infections per person per day", colour = NULL) +
  theme_epiwave
